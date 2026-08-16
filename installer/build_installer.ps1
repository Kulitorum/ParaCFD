# build_installer.ps1 - build the ParaCFD Release GUI (optional) and compile the signed
# Inno Setup installer. Run from anywhere; paths are resolved relative to this script / the repo.
#
#   .\installer\build_installer.ps1                 # compile the SIGNED installer (eToken must be plugged in)
#   .\installer\build_installer.ps1 -Build          # (re)build paracfd-gui Release first, then compile
#   .\installer\build_installer.ps1 -NoSign         # compile an UNSIGNED installer (no token needed)
#   .\installer\build_installer.ps1 -Build -NoSign  # both
#
# The signed build reuses the "SafeNet" named SignTool configured in the Inno IDE on this machine
# (COBOD's Sectigo cert on the SafeNet eToken). A token PIN window pops up during signing; enter it.
# If the token is not present, use -NoSign.
#
# NOTE: kept ASCII-only on purpose so it parses under both Windows PowerShell 5.1 and pwsh 7.
[CmdletBinding()]
param(
    [switch]$Build,   # (re)build paracfd-gui Release before packaging
    [switch]$NoSign   # compile an unsigned installer (passes /DSkipSign to ISCC)
)
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = Split-Path -Parent $here

# Locate ISCC (Inno Setup 6 Command-Line Compiler).
$iscc = @(
    "C:\Program Files (x86)\Inno Setup 6\ISCC.exe",
    "C:\Program Files\Inno Setup 6\ISCC.exe"
) | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $iscc) { throw "ISCC.exe not found. Install Inno Setup 6." }

if ($Build) {
    Write-Host "==> Building paracfd-gui (Release) ..." -ForegroundColor Cyan
    # Free the exe if a previous instance is holding it (it locks during a rebuild).
    taskkill /F /IM paracfd-gui.exe 2>$null | Out-Null
    & cmake --build "$repo\build" --config Release --target paracfd-gui
    if ($LASTEXITCODE -ne 0) { throw "cmake build failed (exit $LASTEXITCODE)" }
}

$rel = "$repo\build\Release"
$exe = "$rel\paracfd-gui.exe"
if (-not (Test-Path $exe)) {
    throw "paracfd-gui.exe not found at $exe. Build it first (pass -Build, or run the CMake Release build)."
}

# Preflight: the runtime closure must be present in build\Release, or the installer ships broken (e.g.
# a missing OCC/Qt DLL = the app won't start on a clean machine). The OCC POST_BUILD deploy only re-runs
# when paracfd-gui RELINKS, so build\Release can silently lack DLLs. If any are missing, -Build forces a
# relink (which redeploys them). NB the 6 OCC 3rdparty DLLs are sourced from the OCC install in the .iss,
# so they are NOT checked here.
$need = @("TKDESTEP.dll","TKService.dll","Qt6Core.dll","Qt6Gui.dll","Qt6Widgets.dll",
    "opengl32sw.dll","platforms\qwindows.dll")
$missing = $need | Where-Object { -not (Test-Path (Join-Path $rel $_)) }
if ($missing) {
    throw ("build\Release is missing runtime files: {0}. Re-run with -Build to force a relink + redeploy." -f ($missing -join ", "))
}

$iss = "$here\ParaCFD-installer.iss"
if ($NoSign) {
    Write-Host "==> Compiling UNSIGNED installer ..." -ForegroundColor Yellow
    & $iscc "/DSkipSign" $iss
}
else {
    Write-Host "==> Compiling SIGNED installer (SafeNet eToken PIN window will appear) ..." -ForegroundColor Cyan
    & $iscc $iss
}
if ($LASTEXITCODE -ne 0) { throw "ISCC failed (exit $LASTEXITCODE)" }

$out = Get-ChildItem "$here\Output\*.exe" -EA SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($out) { Write-Host ("==> Installer written: {0}  ({1} MB)" -f $out.FullName, [math]::Round($out.Length/1MB,1)) -ForegroundColor Green }
else      { Write-Host "==> ISCC reported success but no output .exe found under installer\Output" -ForegroundColor Yellow }
