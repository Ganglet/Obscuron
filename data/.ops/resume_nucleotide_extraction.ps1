# Relaunches the nucleotide panel extraction if it isn't already running.
# Registered as a logon scheduled task (DarkMatterNucleotideResume) so a
# laptop restart doesn't require someone to remember to manually restart
# it -- same pattern as resume_transfer.ps1 for the Week 1 S3 upload.
# Safe to run repeatedly: does nothing if the extraction is already
# running or already finished.

$logFile = "C:\Users\rayya\AppData\Local\Temp\gtdb_nucleotide_extraction.log"
$scriptPath = "E:\dark_matter\data\.ops\run_nucleotide_extraction.sh"
$bash = "F:\Git Files\Git\usr\bin\bash.exe"

$running = Get-CimInstance Win32_Process -Filter "Name='uv.exe' OR Name='python.exe'" |
    Where-Object { $_.CommandLine -like "*extract_panel_nucleotides*" }

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
