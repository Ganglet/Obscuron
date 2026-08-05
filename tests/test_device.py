from darkmatter.device import DevicePlan, detect_device, plan_for_model


def test_detect_device_returns_known_backend():
    device, vram = detect_device()
    assert device in ("cuda", "mps", "cpu")
    if device == "cuda":
        assert vram is not None and vram > 0
    else:
        assert vram is None


def test_plan_quantizes_large_model_on_small_gpu():
    plan = plan_for_model(vram_bf16_gb=9.0, quantize_below_vram_gb=16.0)
    assert isinstance(plan, DevicePlan)
    if plan.device == "cuda" and plan.total_vram_gb is not None and plan.total_vram_gb < 16.0:
        assert plan.quantized is True
        assert plan.quantize_bits == 4


def test_plan_never_quantizes_when_threshold_is_zero():
    plan = plan_for_model(vram_bf16_gb=0.6, quantize_below_vram_gb=0.0)
    assert plan.quantized is False
    assert plan.quantize_bits is None
