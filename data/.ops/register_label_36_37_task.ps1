# Run this from an ELEVATED PowerShell (Run as Administrator).
# Registers the scheduled task that makes the 36->37 boundary scan survive
# a reboot: relaunches the resilient wrapper at logon if the scan isn't
# already running or finished. The wrapper itself already polls every 60s
# and restarts on a stall, so only one task (no separate watchdog) is needed.

$action = New-ScheduledTaskAction -Execute "powershell.exe" -Argument "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File E:\dark_matter\scripts\run_label_36_37_resilient.ps1"
$trigger = New-ScheduledTaskTrigger -AtLogOn
$principal = New-ScheduledTaskPrincipal -UserId "rayya" -LogonType Interactive -RunLevel Limited

Register-ScheduledTask -TaskName "DarkMatterLabel3637Resume" -Action $action -Trigger $trigger -Principal $principal -Description "Resumes the full-panel 36->37 leakage-tightened boundary scan after a reboot/logon if it isn't already running or finished." -Force

Get-ScheduledTask -TaskName "DarkMatterLabel3637Resume" | Select-Object TaskName, State
