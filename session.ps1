#!/usr/bin/env pwsh
<#
.SYNOPSIS
    Start a new Julia development session for the Epicycle workspace.

.DESCRIPTION
    Activates the Epicycle project environment and launches an interactive
    Julia REPL, optionally running setup or tests first.

.PARAMETER Setup
    Run environment setup (Pkg.instantiate + develop all sub-packages) before
    opening the REPL.

.PARAMETER Test
    Run the full local CI test suite instead of opening a REPL.

.PARAMETER Coverage
    Generate coverage reports after testing (implies -Test).

.EXAMPLE
    .\session.ps1
    # Opens a Julia REPL with the Epicycle project activated

.EXAMPLE
    .\session.ps1 -Setup
    # Runs setup_environment.jl then opens the REPL

.EXAMPLE
    .\session.ps1 -Test
    # Runs the local CI test suite
#>

[CmdletBinding()]
param(
    [switch]$Setup,
    [switch]$Test,
    [switch]$Coverage
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Locate the workspace root (the folder containing this script)
# ---------------------------------------------------------------------------
$Root = $PSScriptRoot
if (-not $Root) { $Root = Get-Location }

Write-Host "Epicycle workspace: $Root" -ForegroundColor Cyan

# ---------------------------------------------------------------------------
# Verify julia is on PATH
# ---------------------------------------------------------------------------
if (-not (Get-Command julia -ErrorAction SilentlyContinue)) {
    Write-Error "julia not found on PATH. Please install Julia and ensure it is on your PATH."
    exit 1
}

$JuliaVersion = julia --version 2>&1
Write-Host "Using: $JuliaVersion" -ForegroundColor Green

# ---------------------------------------------------------------------------
# Helper: run a Julia script and stop on failure
# ---------------------------------------------------------------------------
function Invoke-JuliaScript {
    param([string]$Script, [string]$Label)
    Write-Host ""
    Write-Host ("=" * 60) -ForegroundColor DarkGray
    Write-Host "  $Label" -ForegroundColor Yellow
    Write-Host ("=" * 60) -ForegroundColor DarkGray
    & julia --project="$Root" "$Root\$Script"
    if ($LASTEXITCODE -ne 0) {
        Write-Error "$Label failed (exit code $LASTEXITCODE)."
        exit $LASTEXITCODE
    }
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
if ($Coverage) { $Test = $true }

if ($Test) {
    # ---- Full local CI run --------------------------------------------------
    Invoke-JuliaScript "ci\setup_environment.jl"  "Setup Environment"
    Invoke-JuliaScript "ci\build_epicycle.jl"      "Build Epicycle"
    Invoke-JuliaScript "ci\test_epicycle.jl"       "Test Epicycle"
    if ($Coverage) {
        Invoke-JuliaScript "ci\generate_coverage.jl" "Generate Coverage"
    }
    Write-Host ""
    Write-Host "All CI steps completed successfully." -ForegroundColor Green
    exit 0
}

if ($Setup) {
    Invoke-JuliaScript "ci\setup_environment.jl" "Setup Environment"
}

# ---------------------------------------------------------------------------
# Launch interactive REPL with the project activated
# ---------------------------------------------------------------------------
Write-Host ""
Write-Host "Launching Julia REPL (project: $Root)" -ForegroundColor Green
Write-Host 'Type `using Epicycle` to load the package.' -ForegroundColor DarkGray
Write-Host ""

& julia --project="$Root" --banner=yes
