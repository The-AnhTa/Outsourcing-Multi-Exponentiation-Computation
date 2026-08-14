[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(Mandatory = $true)]
    [ValidateSet("run_pippenger", "run_pinkas")]
    [string]$Target,

    [Parameter(Mandatory = $true)]
    [string]$FailureName,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ($RunArguments.Count -ne 4) {
    throw "Usage: --d <dimension> --k <instances>"
}

$values = @{}
for ($i = 0; $i -lt $RunArguments.Count; $i += 2) {
    $key = $RunArguments[$i]
    if ($key -notin @("--d", "-d", "--k", "-k")) {
        throw "Usage: --d <dimension> --k <instances>"
    }
    $normalized = $key.TrimStart("-")
    if ($values.ContainsKey($normalized)) {
        throw "Duplicate argument: $key"
    }
    $values[$normalized] = $RunArguments[$i + 1]
}

$d = 0
$k = 0
if (-not $values.ContainsKey("d") -or
    -not [int]::TryParse($values["d"], [ref]$d) -or
    $d -lt 1 -or $d -gt 30) {
    throw "The dimension must be an integer from 1 through 30."
}
if (-not $values.ContainsKey("k") -or
    -not [int]::TryParse($values["k"], [ref]$k) -or
    $k -lt 1 -or $k -gt 1048576) {
    throw "The instance count must be an integer from 1 through 1048576."
}

$projectRoot = $PSScriptRoot
$buildDirectory = Join-Path $projectRoot "build"
$executable = Join-Path $buildDirectory "$Target.exe"
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
    "cmake --build `"$buildDirectory`" --target $Target -j 4 >nul"

$ErrorActionPreference = "Continue"
$buildOutput = & cmd.exe /d /s /c $nativeCommand 2>&1
$buildExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"
if ($buildExitCode -ne 0) {
    throw "The $Target build failed: $($buildOutput -join [Environment]::NewLine)"
}

$output = & $executable $d $k 2>&1
if ($LASTEXITCODE -ne 0) { throw "$FailureName failed." }
$output
