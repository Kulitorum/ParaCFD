; TEMPLATE for setupvars.iss — machine-specific paths the installer needs.
; Copy this file to `setupvars.iss` (same folder) and adjust the paths for your machine.
; setupvars.iss is committed for this single-dev machine, but keep this template in sync.

; --- The repository root (holds configs\ and Experiments\ shipped as runtime data) --------
#define ProjectDir "C:\CODE\WindCFD"

; --- The CMake Release build output (holds windcfd-gui.exe + the Qt/OpenCascade/3rdparty DLL
;     closure deployed by windeployqt + windcfd_deploy_occ_dlls, and the Qt plugin subdirs) ---
;     Build it first:  cmake --build build --config Release --target windcfd-gui
#define AppBuildDir "C:\CODE\WindCFD\build\Release"

; --- MSVC 2022 (v143) redistributable — installed by [Run]; the app's vcruntime/msvcp are NOT
;     deployed next to the exe, so the target machine needs this. ---------------------------
#define MsvcRedist_Dir "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Redist\MSVC\v143"

; --- ffmpeg.exe bundled next to the app for the video recorder (found via applicationDirPath).
;     Any recent static/essentials ffmpeg build works. Comment out to NOT bundle ffmpeg
;     (recording then needs ffmpeg on PATH or the WINDCFD_FFMPEG env var at run time). ----------
#define FFmpeg_Exe "C:\Users\Micro\Downloads\ffmpeg-2025-12-22-git-c50e5c7778-essentials_build\bin\ffmpeg.exe"

; --- OpenCascade runtime bin dir. The installer sources the OCC 3rdparty DLLs (freetype/brotli/bz2/
;     libpng/zlib) straight from here — they are LOAD-TIME deps of the OCC toolkits, and the build-dir
;     deploy can silently miss them. Adjust to your OpenCascade install (see OCC_BIN_DIR in CMake). ---
#define OCC_BinDir "C:\OpenCASCADE-8.0\build2\win64\vc14\bin"
