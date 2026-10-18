# Runs run_wm_adjusted.R on the lab PC once the composition simulation is
# idle, because a second R/INLA job beside it crashes this PC. Idle means the
# simulation queue mutex is free and no R process from ST 2\tools runs, at
# three checks five minutes apart. Started at logon by the scheduled task
# "SpaPath DLPFC adjusted", so it also resumes after a crash; it exits when
# all 12 x 30 gene sets are finished. The R script resumes from the finished
# gene sets, and a run that ends early (native crash) is retried twice.
$st2 = "C:\Users\lxs888\Documents\ChatGPT\ST 2"
$out = Join-Path $st2 "data\DLPFC_12_sections_20261007\SpaPath_WM_DLPFC30_adjusted_0.7.0"
$script = Join-Path $out "scripts\run_wm_adjusted.R"
$rscript = Join-Path $st2 "tools\R-4.6.1\bin\x64\Rscript.exe"
$log = Join-Path $out "queue.log"
function Note($msg) { Add-Content $log "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $msg" }
function Done-Count { @(Get-ChildItem -Path "$out\1516*\P??" -Directory -ErrorAction SilentlyContinue).Count }

$self = New-Object System.Threading.Mutex($false, "Global\SpaPathDLPFCAdjusted")
if (-not $self.WaitOne(0)) { exit }
if ((Done-Count) -ge 360) { exit }
Note "waiting for the simulation to be idle"
# Give the simulation's own logon task time to start first.
Start-Sleep -Seconds 600

function Is-Idle {
  $queue = New-Object System.Threading.Mutex($false, "Global\SpaPathCompositionQueue")
  $free = $false
  try { $free = $queue.WaitOne(0) } catch [System.Threading.AbandonedMutexException] { $free = $true }
  if ($free) { $queue.ReleaseMutex() }
  $queue.Dispose()
  $r = @(Get-Process -Name Rscript, R, Rterm -ErrorAction SilentlyContinue |
         Where-Object { $_.Path -like "$st2\tools\*" })
  return ($free -and $r.Count -eq 0)
}
$idle = 0
while ($idle -lt 3) {
  if (Is-Idle) { $idle++ } else { $idle = 0 }
  if ($idle -lt 3) { Start-Sleep -Seconds 300 }
}

New-Item -ItemType Directory -Force (Join-Path $out "logs") | Out-Null
for ($attempt = 1; $attempt -le 3 -and (Done-Count) -lt 360; $attempt++) {
  $stamp = Get-Date -Format 'yyyyMMddHHmmss'
  Note "start attempt $attempt ($(Done-Count) of 360 done)"
  $p = Start-Process -FilePath $rscript -ArgumentList "`"$script`"" `
    -RedirectStandardOutput "$out\logs\run_$stamp.out" -RedirectStandardError "$out\logs\run_$stamp.err" `
    -WindowStyle Hidden -PassThru
  $null = $p.Handle
  Start-Sleep -Milliseconds 500
  try { $p.PriorityClass = 'BelowNormal' } catch {}
  $p.WaitForExit()
  Note "end attempt $attempt exit $($p.ExitCode) ($(Done-Count) of 360 done)"
}
$self.ReleaseMutex()
