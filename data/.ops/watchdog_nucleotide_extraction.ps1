# Runs every 10 minutes via scheduled task DarkMatterNucleotideWatchdog.
# Catches the failure mode plain retries can't: the process stays alive
# but gets stuck in a blocked socket read that never times out -- the
# same pattern seen 3x during the Week 1 S3 transfer, right after an SSL
# error, on this exact connection. If the log hasn't grown in 5+ minutes,
# that process is dead weight -- kill it and let the resume logic pick
# back up (extraction is idempotent per-genome, so nothing already
# written is lost, only the in-flight chunk).

$logFile = "C:\Users\rayya\AppData\Local\Temp\gtdb_nucleotide_extraction.log"
$scriptPath = "E:\dark_matter\data\.ops\run_nucleotide_extraction.sh"
$bash = "F:\Git Files\Git\usr\bin\bash.exe"
$staleThresholdMinutes = 5

function Start-Extraction {
    Start-Process -FilePath "cmd.exe" -ArgumentList "/c `"`"$bash`" `"$scriptPath`" >> `"$logFile`" 2>&1`"" -WindowStyle Hidden
}

if (Test-Path $logFile) {
    $tail = Get-Content $logFile -Tail 1 -ErrorAction SilentlyContinue
    if ($tail -like "*ALL DONE*") { exit 0 }
}

$running = Get-CimInstance Win32_Process -Filter "Name='uv.exe' OR Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*extract_panel_nucleotides*" }

if (-not $running) {
    Start-Extraction
    exit 0
}

if (Test-Path $logFile) {
    $staleMinutes = (New-TimeSpan -Start (Get-Item $logFile).LastWriteTime -End (Get-Date)).TotalMinutes
    if ($staleMinutes -gt $staleThresholdMinutes) {
        Get-CimInstance Win32_Process -Filter "Name='python.exe' OR Name='uv.exe' OR Name='bash.exe'" |
            Where-Object { $_.CommandLine -like "*extract_panel_nucleotides*" -or $_.CommandLine -like "*run_nucleotide_extraction*" } |
            ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
        Start-Sleep -Seconds 2
        Start-Extraction
    }
}
