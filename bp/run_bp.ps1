[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

& (Join-Path $PSScriptRoot "run_common.ps1") `
    -Target run_bp `
    -FailureName "BP" `
    -RunArguments $RunArguments
