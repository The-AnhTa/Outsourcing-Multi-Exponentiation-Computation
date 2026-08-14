[CmdletBinding(PositionalBinding = $false)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RunArguments
)

& (Join-Path $PSScriptRoot "run_common.ps1") `
    -Target run_helper_verifier `
    -FailureName "Helper-verifier protocol" `
    -RunArguments $RunArguments
