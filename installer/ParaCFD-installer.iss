; ============================================================================================
; ParaCFD — Inno Setup installer (signed).
;
; Packs the CMake Release build of paracfd-gui.exe together with its full runtime closure:
;   - Qt 6.11 (windeployqt output: Qt6*.dll + the platforms/styles/imageformats/... plugin dirs)
;   - OpenCascade 8.0 STEP toolkits (TK*.dll) + their 3rdparty DLLs (freetype/brotli/bz2/png/zlib)
;   - Qt's software-GL + RHI shader-compiler DLLs (opengl32sw.dll, dxcompiler.dll, dxil.dll)
;   - the app's runtime data (configs\ — the default g1_viewer.json is loaded relative to {app};
;     Experiments\ — the sample STEP models incl. the drop feature's default XStone unit)
;   - a bundled ffmpeg.exe for the MP4 video recorder (found via applicationDirPath at run time)
;   - the MSVC v143 runtime redistributable (run silently post-install)
; CUDA and Jolt are STATICALLY linked (no DLLs to ship); the target machine needs an NVIDIA driver.
;
; SIGNING: reuses the "SafeNet" named SignTool configured in the Inno IDE / registry on this
; machine (COBOD's Sectigo code-signing cert on a SafeNet eToken), identical to cobod-slicer:
;   signtool sign /n "COBOD International A/S" /tr http://timestamp.sectigo.com /td sha256 /fd sha256 $p
; The token must be plugged in (and its password entered when prompted) at COMPILE time. Both the
; installed paracfd-gui.exe (Flags: sign), the uninstaller (SignedUninstaller) and the generated
; setup .exe are signed. If you have no token, comment out SignTool + SignedUninstaller + the
; `sign` flag below to build an UNSIGNED installer.
;
; Machine-specific paths live in setupvars.iss; product identity/version in version.iss.
; Build with build_installer.ps1 (or open this file in the Inno Setup IDE and Compile).
; ============================================================================================

#include "setupvars.iss"
#include "version.iss"

; Pass /DSkipSign to ISCC to build an UNSIGNED installer (e.g. to validate packaging without the
; eToken plugged in). Default (no define) = fully signed. build_installer.ps1 -NoSign sets this.
#ifndef SkipSign
  #define SignFlag "sign"
#else
  #define SignFlag ""
#endif

[Setup]
; A unique AppId for ParaCFD (distinct from cobod-slicer's).
AppId={{CE7E40CF-B660-4C0B-9114-8EF70A9A241C}
AppName={#TargetProduct}
AppVerName={#TargetProduct} v{#VersionNumber}
AppVersion={#VersionNumber}
AppPublisher={#TargetCompany}
DefaultGroupName={#TargetCompany}
DefaultDirName={commonpf}\{#TargetCompany}\{#TargetProduct}
OutputDir=Output
OutputBaseFilename={#TargetName}_v{#VersionNumber}_win{#TargetArch}_installer
UninstallDisplayIcon={app}\{#TargetName}.exe
UninstallDisplayName={#TargetProduct} v{#VersionNumber}
WizardStyle=modern
Compression=lzma2/max
SolidCompression=yes
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
MinVersion=10.0
; Installs into Program Files, so admin rights are required.
PrivilegesRequired=admin

; --- Code signing (SafeNet eToken; see header) -------------------------------------------
; $f expands to the file(s) being signed and is substituted into the named tool's $p parameter.
#ifndef SkipSign
SignTool=SafeNet $f
SignedUninstaller=yes
#endif

[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl"

[CustomMessages]
en.installVcRuntime=Installing the Microsoft Visual C++ runtime ...
en.launchApp=Launch {#TargetProduct}

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

; Wipe the target folder on a re-install so stale DLLs/plugins never linger (must run before [Files]).
[InstallDelete]
Type: filesandordirs; Name: "{app}\*"

[Files]
; --- The main application exe: signed, and listed BEFORE the wildcard + excluded from it so it
;     is packed (and signed) exactly once. -------------------------------------------------
Source: "{#AppBuildDir}\{#TargetName}.exe"; DestDir: "{app}"; Flags: ignoreversion {#SignFlag}

; --- Everything else from the build output: the Qt + OpenCascade + 3rdparty DLLs and the Qt
;     plugin subdirectories (platforms, styles, imageformats, iconengines, networkinformation,
;     tls, generic). Excludes build cruft (pdb/lib/exp/ilk) and the non-shipped test/gate/CLI
;     exes (scour.exe, scour_*_gate.exe, scour_*_tests.exe) + the already-signed app exe. -----
;     The OCC 3rdparty DLLs are excluded here and sourced from the OCC install below (they are LOAD-
;     TIME deps of the OCC toolkits but the build-dir deploy can silently miss them — see below).
Source: "{#AppBuildDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs; \
    Excludes: "*.pdb,*.lib,*.exp,*.ilk,{#TargetName}.exe,scour.exe,scour_*.exe,freetype.dll,brotlicommon.dll,brotlidec.dll,bz2.dll,libpng16.dll,zlib1.dll"

; --- OpenCascade 3rdparty DLLs — sourced STRAIGHT from the OCC install (not the build dir). These are
;     load-time deps of the OCC toolkits (freetype.dll ← TKService.dll); a missing one aborts startup on
;     a clean machine with "freetype.dll not found". No skip flag ⇒ ISCC ERRORS if the OCC bin path is
;     wrong, so we can never again silently ship an installer without them. ----------------------------
Source: "{#OCC_BinDir}\freetype.dll";     DestDir: "{app}"; Flags: ignoreversion
Source: "{#OCC_BinDir}\brotlicommon.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#OCC_BinDir}\brotlidec.dll";    DestDir: "{app}"; Flags: ignoreversion
Source: "{#OCC_BinDir}\bz2.dll";          DestDir: "{app}"; Flags: ignoreversion
Source: "{#OCC_BinDir}\libpng16.dll";     DestDir: "{app}"; Flags: ignoreversion
Source: "{#OCC_BinDir}\zlib1.dll";        DestDir: "{app}"; Flags: ignoreversion

; --- Qt plugin-root pin (see qt.conf). ---------------------------------------------------
Source: "qt.conf"; DestDir: "{app}"; Flags: ignoreversion

; --- Runtime data the app loads by default. configs\g1_viewer.json is the no-argument default
;     (loaded relative to the working dir, which the shortcuts set to {app}). Experiments\ holds
;     the sample STEP models, including the drop feature's default XStone_Decomposed.stp. -----
Source: "{#ProjectDir}\configs\*"; DestDir: "{app}\configs"; Flags: ignoreversion
; Experiments\ is git-ignored (local sample geometry); skip gracefully if it's absent on a machine.
Source: "{#ProjectDir}\Experiments\V000 Code tests\*.stp"; DestDir: "{app}\Experiments\V000 Code tests"; Flags: ignoreversion skipifsourcedoesntexist

; --- Bundled ffmpeg for the MP4 recorder (VideoRecorder::findFfmpeg picks it up from {app}).
;     Comment out FFmpeg_Exe in setupvars.iss to skip bundling; a wrong path skips it too (the
;     recorder then falls back to ffmpeg on PATH or the PARACFD_FFMPEG env var). ------------------
#ifdef FFmpeg_Exe
Source: "{#FFmpeg_Exe}"; DestDir: "{app}"; DestName: "ffmpeg.exe"; Flags: ignoreversion skipifsourcedoesntexist
#endif

; --- MSVC v143 runtime redistributable (deleted after the silent install below). ----------
Source: "{#MsvcRedist_Dir}\vc_redist.{#TargetArch}.exe"; DestDir: "{app}"; Flags: deleteafterinstall

[Icons]
Name: "{group}\{#TargetProduct}"; Filename: "{app}\{#TargetName}.exe"; WorkingDir: "{app}"
Name: "{group}\Uninstall {#TargetProduct}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\{#TargetProduct}"; Filename: "{app}\{#TargetName}.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Run]
; Install the VC++ runtime silently, then optionally launch the app.
Filename: "{app}\vc_redist.{#TargetArch}.exe"; Parameters: "/install /quiet /norestart"; StatusMsg: "{cm:installVcRuntime}"
Filename: "{app}\{#TargetName}.exe"; Description: "{cm:launchApp}"; WorkingDir: "{app}"; Flags: nowait postinstall skipifsilent
