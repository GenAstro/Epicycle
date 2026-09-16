param([int] $Runs = 3)
$ErrorActionPreference = "Continue"

Write-Host "`n############ Scenario A: Pkg.test() (fresh Julia each time) ############"
$pkgTest = @()
for ($i = 1; $i -le $Runs; $i++) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $out = & julia --project=. -e "using Pkg; Pkg.test()" 2>&1
    $sw.Stop()
    $summary = ($out | Select-String -Pattern "Package\s+\|" | Select-Object -First 1)
    $summary = if ($summary) { $summary.ToString().Trim() } else { "(no summary line)" }
    $pkgTest += [pscustomobject]@{
        Run          = $i
        WallClockSec = [math]::Round($sw.Elapsed.TotalSeconds, 2)
        Summary      = $summary
    }
    Write-Host ("  run {0}: {1}s  ({2})" -f $i, [math]::Round($sw.Elapsed.TotalSeconds,2), $summary)
}

Write-Host "`n############ Scenario B: bench_tests.jl (fresh Julia, single run) ############"
$coldSw = [Diagnostics.Stopwatch]::StartNew()
$coldOut = & julia --project=./benchenv .\bench_tests.jl 2>&1
$coldSw.Stop()
$coldOut | ForEach-Object { Write-Host $_ }
Write-Host ("  wall-clock: {0}s" -f [math]::Round($coldSw.Elapsed.TotalSeconds,2))

Write-Host "`n############ Scenario C: warm loop (one Julia session, N includes) ############"
$expr = ''
for ($i = 1; $i -le $Runs; $i++) { $expr += 'include("bench_tests.jl"); ' }
$warmSw = [Diagnostics.Stopwatch]::StartNew()
$warmOut = & julia --project=./benchenv -e $expr 2>&1
$warmSw.Stop()
$warmOut | ForEach-Object { Write-Host $_ }
Write-Host ("  total wall-clock (session with $Runs includes): {0}s" -f [math]::Round($warmSw.Elapsed.TotalSeconds,2))

Write-Host "`n===================== Pkg.test summary ====================="
$pkgTest | Format-Table -AutoSize
