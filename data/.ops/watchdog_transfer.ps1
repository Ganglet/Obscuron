# Runs every 10 minutes via scheduled task DarkMatterGtdbWatchdog.
# Catches the failure mode plain retries can't: the transfer process stays
# alive but gets stuck in a blocked socket read that never times out
# (seen 3x on this connection, always right after an SSL error - looks
# like local HTTPS-inspection security software occasionally corrupting a
# stream badly enough that neither requests' nor boto3's timeouts fire).
# If the log hasn't grown in 5+ minutes, that process is dead weight -
# kill it and let the normal resume logic pick back up from S3's state.

$logFile = "C:\Users\rayya\AppData\Local\Temp\gtdb_transfer.log"
$scriptPath = "E:\dark_matter\data\.ops\run_gtdb_transfer.sh"
$bash = "F:\Git Files\Git\usr\bin\bash.exe"
$staleThresholdMinutes = 5

function Start-Transfer {
    Start-Process -FilePath "cmd.exe" -ArgumentList "/c `"`"$bash`" `"$scriptPath`" >> `"$logFile`" 2>&1`"" -WindowStyle Hidden
}

if (Test-Path $logFile) {
    $tail = Get-Content $logFile -Tail 1 -ErrorAction SilentlyContinue
    if ($tail -like "*ALL DONE*") { exit 0 }
}

$running = Get-CimInstance Win32_Process -Filter "Name='uv.exe' OR Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*stream_gtdb_to_s3*" }

if (-not $running) {
    Start-Transfer
    exit 0
}

if (Test-Path $logFile) {
    $staleMinutes = (New-TimeSpan -Start (Get-Item $logFile).LastWriteTime -End (Get-Date)).TotalMinutes
    if ($staleMinutes -gt $staleThresholdMinutes) {
        Get-CimInstance Win32_Process -Filter "Name='python.exe' OR Name='uv.exe' OR Name='bash.exe'" |
            Where-Object { $_.CommandLine -like "*stream_gtdb*" -or $_.CommandLine -like "*run_gtdb_transfer*" } |
            ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
        Start-Sleep -Seconds 2
        Start-Transfer
    }
}
