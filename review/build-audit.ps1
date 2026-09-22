$ErrorActionPreference = 'Stop'
Set-Location (Split-Path $PSScriptRoot -Parent)
$auditVs = & 'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe' -version '[17.0,18.0)' -property installationPath
if (-not $auditVs) { throw 'Visual Studio 2022 is required for this local diagnostic build.' }
$auditCuda = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$auditOcc = 'C:\OpenCASCADE-8.0\build2\win64\vc14\lib'
$auditOutput = 'build-paraglider-ui\Release'
$auditToolkits = 'TKernel','TKMath','TKG2d','TKG3d','TKGeomBase','TKGeomAlgo','TKBRep','TKTopAlgo','TKMesh','TKPrim','TKBO','TKBool','TKShHealing','TKXSBase','TKDESTEP'
$auditLibraries = ($auditToolkits | ForEach-Object { '"' + $auditOcc + '\' + $_ + '.lib"' }) -join ' '
$auditCommand = 'call "' + $auditVs + '\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /std:c++20 /EHsc /MD /O2 /I src /I src\3rdparty /I "' + $auditCuda + '\include" review\physics_audit_probe.cpp /Fo"' + $auditOutput + '\physics_audit_probe.obj" /Fe"' + $auditOutput + '\physics_audit_probe.exe" /link /NODEFAULTLIB:LIBCMT "' + $auditOutput + '\libparacfd.lib" "' + $auditOutput + '\paracfd_geometry.lib" "' + $auditCuda + '\lib\x64\cudart_static.lib" ' + $auditLibraries
& cmd.exe /d /c $auditCommand
if ($LASTEXITCODE -ne 0) { throw 'Diagnostic compilation failed.' }
