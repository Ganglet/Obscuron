# Relaunches the GTDB->S3 transfer if it isn't already running.
# Registered as a logon scheduled task (DarkMatterGtdbResume) so a laptop
# restart doesn't require someone to remember to manually restart it.
# Safe to run repeatedly: does nothing if the transfer is already running
# or already finished.

$logFile = "C:\Users\rayya\AppData\Local\Temp\gtdb_transfer.log"
$scriptPath = "E:\dark_matter\data\.ops\run_gtdb_transfer.sh"
$bash = "F:\Git Files\Git\usr\bin\bash.exe"

$running = Get-CimInstance Win32_Process -Filter "Name='uv.exe' OR Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*stream_gtdb_to_s3*" }

if ($running) {
    exit 0
}

if (Test-Path $logFile) {
    $tail = Get-Content $logFile -Tail 1 -ErrorAction SilentlyContinue
    if ($tail -like "*ALL DONE*") {
        exit 0
    }
}

Start-Process -FilePath "cmd.exe" -ArgumentList "/c `"`"$bash`" `"$scriptPath`" >> `"$logFile`" 2>&1`"" -WindowStyle Hidden
