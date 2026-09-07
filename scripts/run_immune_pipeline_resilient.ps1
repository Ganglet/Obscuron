# Resilient driver for the Phase-3 Track-2 immune pipeline (reference re-embed
# -> dark-query re-embed -> sweep). Survives laptop restarts (run via a
# Task Scheduler ONLOGON trigger, not tied to any terminal/session), WiFi
# drops (offline mode -- both models are already cached locally), and a
# hung/stuck step (a watchdog restarts it if its checkpoint stops advancing).
# Each underlying python step is independently checkpointed (P1-D10 pattern),
# so a kill-and-retry loses at most one chunk, not the whole run.
#
# Manual re-run is always safe -- every step no-ops if its output already
# exists.

$ErrorActionPreference = 'Continue'
Set-Location "E:\dark_matter"

$env:HF_HOME = "E:\dark_matter\.cache\huggingface"
$env:HF_HUB_OFFLINE = "1"
$env:TRANSFORMERS_OFFLINE = "1"

$Proc  = "E:\dark_matter\data\processed\gtdb_R207"
$Results = "E:\dark_matter\results"
$LogFile = "$Proc\pipeline_run.log"
$DoneMarker = "$Proc\pipeline_DONE"
$StaleMinutes = 25
$PollSeconds = 60
$MaxRetries = 30

function Log($msg) {
    $line = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  $msg"
    Write-Output $line
    Add-Content -Path $LogFile -Value $line
}

function Run-Step($name, $scriptPath, $scriptArgs, $outputCheck, $checkpointCheck) {
    if (Test-Path $outputCheck) {
        Log "[$name] already done ($outputCheck exists) -- skipping"
        return $true
    }
    $stdOut = "$Proc\$name.stdout.log"
    for ($attempt = 1; $attempt -le $MaxRetries; $attempt++) {
        Log "[$name] attempt $attempt/$MaxRetries -- launching"
        Remove-Item $stdOut -ErrorAction SilentlyContinue

        # Start-Job (out-of-process) instead of Start-Process: Start-Process
        # -NoNewWindow combined with dual stdout/stderr redirection is an
        # unreliable combo on Windows PowerShell (HasExited can go true
        # within seconds of a healthy long-running child) -- Start-Job's
        # own stream capture does not have that failure mode.
        $job = Start-Job -ScriptBlock {
            param($work, $script, $args_)
            Set-Location $work
            & uv run python $script @args_ *>&1
        } -ArgumentList "E:\dark_matter", $scriptPath, $scriptArgs

        $lastMtime = Get-Date
        $lastSize = -1
        while ($job.State -eq 'Running') {
            Start-Sleep -Seconds $PollSeconds
            Receive-Job -Job $job -Keep 2>&1 | Out-File -FilePath $stdOut -Encoding utf8
            if (Test-Path $checkpointCheck) {
                $size = (Get-Item $checkpointCheck).Length
            } else {
                $size = (Get-Item $stdOut -ErrorAction SilentlyContinue).Length
            }
            if ($size -ne $lastSize) { $lastSize = $size; $lastMtime = Get-Date }
            $staleFor = (Get-Date) - $lastMtime
            if ($staleFor.TotalMinutes -ge $StaleMinutes) {
                Log "[$name] no progress for $StaleMinutes min -- killing and retrying"
                Stop-Job -Job $job
                break
            }
        }
        Receive-Job -Job $job -Keep 2>&1 | Out-File -FilePath $stdOut -Encoding utf8
        $state = $job.State
        Remove-Job -Job $job -Force -ErrorAction SilentlyContinue

        if (Test-Path $outputCheck) {
            Log "[$name] finished cleanly (job state=$state, output present)"
            return $true
        }
        Log "[$name] attempt $attempt ended (job state=$state, output not yet present) -- backing off 30s before retry"
        Start-Sleep -Seconds 30
    }
    Log "[$name] FAILED after $MaxRetries attempts -- giving up, check $stdOut"
    return $false
}

if (Test-Path $DoneMarker) {
    Log "pipeline already complete ($DoneMarker exists) -- nothing to do"
    exit 0
}

New-Item -ItemType Directory -Force -Path $Results | Out-Null
Log "=== resilient immune pipeline run starting ==="

$ok1 = Run-Step "reembed_reference" "scripts/reembed_reference_immune.py" `
    @("--layers", "33,22", "--fp32", "--batch-size", "16", "--chunk-size", "1000") `
    "$Proc\esm2_reembed_L33_22.npz" "$Proc\esm2_reembed_L33_22.ckpt.npz"
if (-not $ok1) { Log "aborting: reference re-embed did not complete"; exit 1 }

$ok2 = Run-Step "reembed_queries" "scripts/reembed_dark_queries.py" `
    @("--layers", "33,22", "--fp32", "--batch-size", "16", "--chunk-size", "1000", "--n-dark-negative", "10000") `
    "$Proc\esm2_reembed_queries_L33_22.npz" "$Proc\esm2_reembed_queries_L33_22.ckpt.npz"
if (-not $ok2) { Log "aborting: dark-query re-embed did not complete"; exit 1 }

$ok3 = Run-Step "immune_sweep" "scripts/immune_sweep.py" @() `
    "$Results\immune_scale_summary.json" "$Results\immune_scale_summary.json"
if (-not $ok3) { Log "aborting: sweep did not complete"; exit 1 }

New-Item -ItemType File -Force -Path $DoneMarker | Out-Null
Log "=== pipeline complete ==="
