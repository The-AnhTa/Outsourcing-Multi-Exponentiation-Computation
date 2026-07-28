[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ($RunArguments.Count -ne 2 -or
    $RunArguments[0] -notin @("--d", "-d")) {
    throw "Usage: .\run_helper_verifier.ps1 --d <dimension>"
}

$d = 0
if (-not [int]::TryParse($RunArguments[1], [ref]$d) -or
    $d -lt 1 -or $d -gt 20) {
    throw "The dimension must be an integer from 1 through 20."
}

$projectRoot = $PSScriptRoot
$buildDirectory = Join-Path $projectRoot "build"
$executable = Join-Path $buildDirectory "run_helper_verifier.exe"
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path -LiteralPath $vswhere)) { throw "vswhere.exe was not found." }
$installation = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $installation) { throw "Visual Studio C++ Build Tools were not found." }
$vcvars = Join-Path $installation "VC\Auxiliary\Build\vcvars64.bat"

$configureCommand = ""
if (-not (Test-Path -LiteralPath (Join-Path $buildDirectory "build.ninja"))) {
    $configureCommand =
        "cmake -S `"$projectRoot`" -B `"$buildDirectory`" -G Ninja " +
        "-DCMAKE_BUILD_TYPE=Release && "
}

$nativeCommand =
    "call `"$vcvars`" >nul && " +
    $configureCommand +
    "cmake --build `"$buildDirectory`" --target run_helper_verifier -j 4 >nul"

$ErrorActionPreference = "Continue"
$buildOutput = & cmd.exe /d /s /c $nativeCommand 2>&1
$buildExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($buildExitCode -ne 0) {
    throw "The helper-verifier build failed: $($buildOutput -join [Environment]::NewLine)"
}

$output = & $executable $d 2>&1
if ($LASTEXITCODE -ne 0) {
    throw "The helper-verifier protocol failed."
}

$output
