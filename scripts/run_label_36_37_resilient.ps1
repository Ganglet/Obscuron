# Resilient driver for the full-panel 36->37 leakage-tightened boundary scan
# (Layer-4 Genos-m convergence-axis follow-up). Survives laptop restarts (run
# via a Task Scheduler ONLOGON trigger, not tied to any terminal/session) and
# a hung/stuck step (a watchdog restarts it if the output CSV stops growing).
# No network/wifi dependency at all -- this is pure local hmmscan against
# already-downloaded Pfam HMM files, so a wifi drop cannot affect it; only a
# process kill (shutdown, crash) matters, and the underlying python script is
# already checkpointed at batch granularity (P1-D10 pattern) -- a kill-and-
# retry loses at most one ~13-minute batch, not the whole ~5-hour run.
#
# Manual re-run is always safe -- the underlying script skips every genome
# already present in the output CSV.

$ErrorActionPreference = 'Continue'
Set-Location "E:\dark_matter"

$Proc = "E:\dark_matter\data\processed\gtdb_R207"
$LogFile = "$Proc\label_36_37_pipeline.log"
$OutputCsv = "$Proc\panel_protein_labels_36_37.csv"
$DoneMarker = "$Proc\panel_protein_labels_36_37.DONE"
$StaleMinutes = 20
$PollSeconds = 60
$MaxRetries = 30

function Log($msg) {
    $line = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  $msg"
    Write-Output $line
    Add-Content -Path $LogFile -Value $line
}

if (Test-Path $DoneMarker) {
    Log "scan already complete ($DoneMarker exists) -- nothing to do"
    exit 0
}

Log "=== resilient 36->37 boundary scan starting ==="

for ($attempt = 1; $attempt -le $MaxRetries; $attempt++) {
    Log "attempt $attempt/$MaxRetries -- launching"
    $stdOut = "$Proc\label_36_37.stdout.log"
    Remove-Item $stdOut -ErrorAction SilentlyContinue

    $job = Start-Job -ScriptBlock {
        param($work)
        Set-Location $work
        & uv run python scripts/label_panel_proteins_36_37.py --batch-size 25 --cpus 16 *>&1
    } -ArgumentList "E:\dark_matter"

    $lastMtime = Get-Date
    $lastSize = -1
    while ($job.State -eq 'Running') {
        Start-Sleep -Seconds $PollSeconds
        Receive-Job -Job $job -Keep 2>&1 | Out-File -FilePath $stdOut -Encoding utf8
        if (Test-Path $OutputCsv) {
            $size = (Get-Item $OutputCsv).Length
        } else {
            $size = (Get-Item $stdOut -ErrorAction SilentlyContinue).Length
        }
        if ($size -ne $lastSize) { $lastSize = $size; $lastMtime = Get-Date }
        $staleFor = (Get-Date) - $lastMtime
        if ($staleFor.TotalMinutes -ge $StaleMinutes) {
            Log "no progress for $StaleMinutes min -- killing and retrying"
            Stop-Job -Job $job
            break
        }
    }
    Receive-Job -Job $job -Keep 2>&1 | Out-File -FilePath $stdOut -Encoding utf8
    $state = $job.State
    Remove-Job -Job $job -Force -ErrorAction SilentlyContinue

    if (Test-Path $DoneMarker) {
        Log "finished cleanly (job state=$state)"
        exit 0
    }
    Log "attempt $attempt ended (job state=$state, not yet done) -- backing off 30s before retry"
    Start-Sleep -Seconds 30
}

Log "FAILED after $MaxRetries attempts -- giving up, check $stdOut"
exit 1
