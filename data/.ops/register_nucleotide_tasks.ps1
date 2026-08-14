# Run this from an ELEVATED PowerShell (Run as Administrator).
# Registers the two scheduled tasks that make nucleotide extraction
# survive a reboot/wifi drop: one relaunches it at logon if it isn't
# already running or finished, the other checks every 10 minutes for a
# hung (not crashed) process and restarts it.

$resumeAction = New-ScheduledTaskAction -Execute "powershell.exe" -Argument "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File E:\dark_matter\data\.ops\resume_nucleotide_extraction.ps1"
$resumeTrigger = New-ScheduledTaskTrigger -AtLogOn
$principal = New-ScheduledTaskPrincipal -UserId "rayya" -LogonType Interactive -RunLevel Limited

Register-ScheduledTask -TaskName "DarkMatterNucleotideResume" -Action $resumeAction -Trigger $resumeTrigger -Principal $principal -Description "Resumes the GTDB nucleotide panel extraction after a reboot/logon if it isn't already running or finished." -Force

$watchdogAction = New-ScheduledTaskAction -Execute "powershell.exe" -Argument "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File E:\dark_matter\data\.ops\watchdog_nucleotide_extraction.ps1"
$watchdogTrigger = New-ScheduledTaskTrigger -Once -At (Get-Date) -RepetitionInterval (New-TimeSpan -Minutes 10) -RepetitionDuration (New-TimeSpan -Days 3650)

Register-ScheduledTask -TaskName "DarkMatterNucleotideWatchdog" -Action $watchdogAction -Trigger $watchdogTrigger -Principal $principal -Description "Every 10 min: kills and restarts the nucleotide extraction if its log has gone stale (stuck socket read), or starts it if not running." -Force

Get-ScheduledTask -TaskName "DarkMatterNucleotide*" | Select-Object TaskName, State
