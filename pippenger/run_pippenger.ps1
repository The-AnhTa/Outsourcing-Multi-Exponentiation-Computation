[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

& (Join-Path $PSScriptRoot "run_common.ps1") `
    -Target run_pippenger `
    -FailureName "Multi-instance Pippenger" `
    -RunArguments $RunArguments
