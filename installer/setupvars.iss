; setupvars.iss — machine-specific paths (MH's workstation). See TEMPLATE_setupvars.iss.
; Committed because this is a single-dev machine; edit both together if a path moves.

; Repository root (ships configs\ + Experiments\ as runtime data).
#define ProjectDir "C:\CODE\ParaCFD"

; CMake Release build output (paracfd-gui.exe + the full deployed DLL/plugin closure).
; Build first:  cmake --build build --config Release --target paracfd-gui
#define AppBuildDir "C:\CODE\ParaCFD\build\Release"

; MSVC 2022 (v143) redistributable, run silently post-install.
#define MsvcRedist_Dir "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Redist\MSVC\v143"

; ffmpeg.exe bundled for the video recorder (VideoRecorder::findFfmpeg picks it up from {app}).
#define FFmpeg_Exe "C:\Users\Micro\Downloads\ffmpeg-2025-12-22-git-c50e5c7778-essentials_build\bin\ffmpeg.exe"

; OpenCascade runtime bin dir. The installer sources the OCC 3rdparty DLLs (freetype/brotli/bz2/
; libpng/zlib) straight from here rather than trusting the build-dir POST_BUILD deploy (which only
; re-runs when paracfd-gui relinks, so build/Release can silently lack them — freetype.dll is a
; LOAD-TIME dep of TKService.dll, so a missing one = the app won't even start on a clean machine).
#define OCC_BinDir "C:\OpenCASCADE-8.0\build2\win64\vc14\bin"
