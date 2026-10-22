param()
$ErrorActionPreference = "Continue"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

$ruffOut = uv run ruff check src tests 2>&1
$ruffOk = $LASTEXITCODE -eq 0

$mypyOut = uv run mypy src 2>&1
$mypyOk = $LASTEXITCODE -eq 0

$pytestOut = uv run pytest -q 2>&1
$pytestOk = $LASTEXITCODE -eq 0

$status = "STATUS ruff=$(if ($ruffOk) {'OK'} else {'FAIL'}) mypy=$(if ($mypyOk) {'OK'} else {'FAIL'}) pytest=$(if ($pytestOk) {'OK'} else {'FAIL'})"

$log = @()
$log += $status
$log += "---- ruff ----"
$log += $ruffOut
$log += "---- mypy ----"
$log += $mypyOut
$log += "---- pytest ----"
$log += $pytestOut
$log -join "`n" | Out-File -FilePath "$root\status_log.txt" -Encoding utf8

Write-Output $status
if (-not ($ruffOk -and $mypyOk -and $pytestOk)) { exit 1 }
