[CmdletBinding()]
param(
    [ValidateRange(0, 3600)]
    [int]$DelaySeconds = 0
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$runner = Join-Path $PSScriptRoot "run_pinkas.ps1"
$dimensions = @(8, 9, 10, 11, 12)
$instanceCounts = @(1, 2, 4, 8)

if (-not (Test-Path -LiteralPath $runner)) {
    throw "Pinkas runner was not found: $runner"
}

foreach ($d in $dimensions) {
    foreach ($k in $instanceCounts) {
        Write-Output "d=$d k=$k"
        & $runner --d $d --k $k
        if ($DelaySeconds -gt 0) {
            Start-Sleep -Seconds $DelaySeconds
        }
    }
}
