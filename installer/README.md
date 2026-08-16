# ParaCFD installer (Inno Setup, signed)

Builds a signed Windows installer that packs `paracfd-gui.exe` with its **entire** runtime closure
so it runs on a clean machine (no Qt / OpenCascade / dev tools required — only an NVIDIA driver).

## What it packs
- **App**: `paracfd-gui.exe` (signed).
- **Qt 6.11**: `Qt6*.dll` + the plugin subdirs (`platforms/`, `styles/`, `imageformats/`,
  `iconengines/`, `networkinformation/`, `tls/`, `generic/`) as deployed by `windeployqt`, plus
  `opengl32sw.dll`, `dxcompiler.dll`, `dxil.dll` and a `qt.conf` pinning the plugin root.
- **OpenCascade 8.0** STEP toolkits (`TK*.dll`) and their 3rdparty DLLs
  (`freetype`, `brotli*`, `bz2`, `libpng16`, `zlib1`) — deployed by `paracfd_deploy_occ_dlls`.
- **Runtime data**: `configs\` (the no-argument default `g1_viewer.json` is loaded relative to
  `{app}`, which the shortcuts set as the working dir) and `Experiments\V000 Code tests\*.stp`
  (the sample STEP models, incl. the drop feature's default `XStone_Decomposed.stp`).
- **ffmpeg.exe** for the MP4 video recorder — bundled into `{app}` and found via
  `VideoRecorder::findFfmpeg()` (`applicationDirPath`), so recording works out of the box.
- **MSVC v143 runtime** (`vc_redist.x64.exe`), run silently post-install.

CUDA 13.1 and Jolt are **statically linked** — no DLLs to ship. The target machine only needs a
compatible **NVIDIA GPU driver** (RTX / CUDA-capable). The build is **pre-calibration / qualitative**
(see `CLAUDE.md`, `M6_REPORT.md`).

## Signing
Reuses the **`SafeNet`** named SignTool configured in the Inno Setup IDE / registry on this machine —
the same COBOD Sectigo code-signing certificate on a SafeNet eToken that cobod-slicer uses:

```
signtool sign /n "COBOD International A/S" /tr http://timestamp.sectigo.com /td sha256 /fd sha256 $p
```

The **eToken must be plugged in** (and its password entered when prompted) at *compile* time. The
installed `paracfd-gui.exe`, the uninstaller (`SignedUninstaller`) and the generated setup `.exe` are
all signed. To register the tool on a new machine: Inno IDE → *Tools → Configure Sign Tools…* → add
one named `SafeNet` with the command above.

## Build it
```powershell
# 1) Build the Release GUI (once; or pass -Build to the script below)
cmake --build build --config Release --target paracfd-gui

# 2) Compile the installer
.\installer\build_installer.ps1            # SIGNED (eToken plugged in)  → installer\Output\*.exe
.\installer\build_installer.ps1 -Build     # rebuild paracfd-gui first, then compile
.\installer\build_installer.ps1 -NoSign    # UNSIGNED (validate packaging without the token)
```
Or open `ParaCFD-installer.iss` in the Inno Setup IDE and press *Compile*.

## Files
- `ParaCFD-installer.iss` — the installer script.
- `version.iss` — product name/company/version (bump `VersionNumber` per release).
- `setupvars.iss` — machine-specific paths (build dir, MSVC redist, ffmpeg). Edit for your box;
  `TEMPLATE_setupvars.iss` is the documented template.
- `qt.conf` — pins Qt's plugin root to `{app}`.
- `build_installer.ps1` — build + compile helper.
- `Output\` — where the compiled installer lands (git-ignored).
