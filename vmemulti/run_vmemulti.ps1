[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ($RunArguments.Count -ne 4) {
    throw "Usage: .\run_vmemulti.ps1 --d <dimension> --m <instances>"
}

$values = @{}
for ($i = 0; $i -lt $RunArguments.Count; $i += 2) {
    $key = $RunArguments[$i]
    if ($key -notin @("--d", "-d", "--m", "-m") -or
        $values.ContainsKey($key.TrimStart("-"))) {
        throw "Usage: .\run_vmemulti.ps1 --d <dimension> --m <instances>"
    }
    $values[$key.TrimStart("-")] = $RunArguments[$i + 1]
}

$d = 0
$m = 0
if (-not $values.ContainsKey("d") -or
    -not [int]::TryParse($values["d"], [ref]$d) -or
    $d -lt 1 -or $d -gt 12) {
    throw "The dimension must be an integer from 1 through 12."
}
if (-not $values.ContainsKey("m") -or
    -not [int]::TryParse($values["m"], [ref]$m) -or
    $m -lt 1 -or $m -gt 1024) {
    throw "The instance count must be an integer from 1 through 1024."
}

$projectRoot = $PSScriptRoot
$buildDirectory = Join-Path $projectRoot "build"
$executable = Join-Path $buildDirectory "run_vmemulti.exe"
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
    "cmake --build `"$buildDirectory`" --target run_vmemulti -j 4 >nul"

$ErrorActionPreference = "Continue"
$buildOutput = & cmd.exe /d /s /c $nativeCommand 2>&1
$buildExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($buildExitCode -ne 0) {
    throw "The run_vmemulti build failed: $($buildOutput -join [Environment]::NewLine)"
}

$output = & $executable $d $m 2>&1
if ($LASTEXITCODE -ne 0) {
    throw "Multi-instance VME failed."
}

$output
