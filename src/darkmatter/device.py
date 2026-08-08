"""Hardware detection and dtype/quantization policy.

Centralizes the "which device, which precision, quantize or not" decision so
every embedding backend applies the same VRAM-aware policy instead of each
guessing independently.
"""

from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass(frozen=True)
class DevicePlan:
    device: str          # "cuda", "mps", or "cpu"
    total_vram_gb: float | None   # None for mps/cpu (not queryable the same way)
    quantize_bits: int | None     # None, 4, or 8
    torch_dtype: torch.dtype

    @property
    def quantized(self) -> bool:
        return self.quantize_bits is not None


def detect_device() -> tuple[str, float | None]:
    if torch.cuda.is_available():
        total_bytes = torch.cuda.get_device_properties(0).total_memory
        return "cuda", total_bytes / (1024**3)
    if torch.backends.mps.is_available():
        return "mps", None
    return "cpu", None


def plan_for_model(vram_bf16_gb: float, quantize_below_vram_gb: float) -> DevicePlan:
    """Decide device/dtype/quantization for a model given its bf16 footprint.

    quantize_below_vram_gb=0 means "never quantize this model" (e.g. it's
    already small enough that quantization isn't worth the accuracy cost).

    Quantizes to 4-bit (not 8-bit) when needed: MoE models keep every expert
    resident in memory regardless of how few activate per token, so on an
    8GB card even 8-bit can be too tight once CUDA context and activations
    are accounted for — confirmed empirically on Genos-m (4.7B total
    params), where 8-bit required CPU offload that then failed on this
    machine's constrained pagefile. 4-bit (~2.5GB for Genos-m) fits with
    headroom to spare.
    """
    device, total_vram_gb = detect_device()

    if device == "cuda":
        should_quantize = (
            quantize_below_vram_gb > 0 and total_vram_gb is not None
            and total_vram_gb < quantize_below_vram_gb
        )
        return DevicePlan(
            device=device,
            total_vram_gb=total_vram_gb,
            quantize_bits=4 if should_quantize else None,
            torch_dtype=torch.bfloat16,
        )

    if device == "mps":
        # bitsandbytes doesn't support MPS; fall back to fp16, no quantization.
        return DevicePlan(device=device, total_vram_gb=None, quantize_bits=None, torch_dtype=torch.float16)

    return DevicePlan(device=device, total_vram_gb=None, quantize_bits=None, torch_dtype=torch.float32)


def peak_vram_gb() -> float | None:
    if not torch.cuda.is_available():
        return None
    return torch.cuda.max_memory_allocated() / (1024**3)


def from_pretrained_kwargs(plan: DevicePlan) -> dict:
    """`transformers.AutoModel.from_pretrained` kwargs implementing a DevicePlan.

    transformers>=5 removed the old `load_in_8bit=True` shorthand in favor of
    an explicit BitsAndBytesConfig — centralized here so both embedding
    backends get quantized loading the same (correct) way.
    """
    kwargs: dict = {"dtype": plan.torch_dtype}
    if plan.quantized:
        from transformers import BitsAndBytesConfig

        if plan.quantize_bits == 4:
            kwargs["quantization_config"] = BitsAndBytesConfig(
                load_in_4bit=True,
                bnb_4bit_quant_type="nf4",
                bnb_4bit_compute_dtype=plan.torch_dtype,
            )
        else:
            kwargs["quantization_config"] = BitsAndBytesConfig(
                load_in_8bit=True, llm_int8_enable_fp32_cpu_offload=True
            )
        # "auto" device_map's memory-estimation heuristic proved overly
        # conservative here (reported overflow on a model that measures
        # well under budget once quantized) — pin everything to the single
        # GPU explicitly instead of trusting the auto-dispatch estimate.
        kwargs["device_map"] = {"": 0}
    return kwargs


def inference_device(plan: DevicePlan, model) -> str:
    """Device to move input tensors to. With device_map="auto" (quantized
    loading) the model's own placement is authoritative, not plan.device."""
    return model.device if plan.quantized else plan.device
