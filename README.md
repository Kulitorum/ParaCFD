# ParaCFD

GPU-accelerated experimental CFD for paraglider aerodynamics, written in C++20 and CUDA with a Qt 6 / OpenGL desktop interface and Open CASCADE STEP import.

> **Current status: the app builds and runs, but it does not reliably calculate the lift of the supplied paraglider. That is the main purpose of the project, and it remains broken.** Force plots, glide estimates, and successful individual diagnostics should not be treated as validated aerodynamic predictions.

This project is being opened up so others can investigate the numerical problems and take the solver further. There is a working interface, CAD import, a GPU solver, example geometry, and a substantial record of attempted fixes. There is still no accepted end-to-end solution for the main wing case. Contributions from people with experience in CFD, cut-cell methods, pressure–velocity coupling, and lifting flows are especially welcome.

## What is here

- STEP import and tessellation through Open CASCADE, with a closed-solid external-flow model.
- Static Cartesian adaptive mesh refinement (AMR), embedded boundaries, and cut-cell geometry.
- CUDA momentum transport, pressure projection, and surface-load integration.
- An interactive Qt viewer with velocity/pressure slices, surface pressure, flow arrows, tracers, and recorded playback.
- Experimental angle sweeps and glide/trim calculations in the UI.
- PlanB paraglider geometry, NACA airfoils, smaller geometry fixtures, and command-line investigation tools.

These are implemented capabilities, not claims of physical correctness. The current model uses fixed, rigid geometry; it does not simulate canopy inflation or coupled fabric deformation. CUDA is required: there is no supported CPU-only solver build.

## The problem to solve

The main case is the supplied PlanB wing at its existing neutral line trim. The goal is credible lift in that configuration, consistent with the intended flight condition. Changing its angle until a force plot looks plausible does not resolve the underlying problem.

Earlier investigations found problems in geometry classification, pressure/momentum consistency, surface-force integration, and small-cell stabilization. Several local defects have been repaired, but those repairs have not established correct whole-wing lift. Some simpler airfoil comparisons and operator checks have passed historically; they do not validate the main case.

Start with these documents, in this order:

1. [Latest recovery attempt](review/FINAL_ATTEMPT.md) — the most recent recorded changes, results, and unresolved lift problem.
2. [Simulation review](SIMULATION_REVIEW.md) — the earlier numerical audit and its reasoning; some findings were addressed by subsequent work.
3. [Completion plan](COMPLETION_PLAN.md) — proposed next steps and the intended physical acceptance case.
4. [Review directory](review/README.md) — diagnostic scripts, intermediate reports, and saved measurements.

Reports are a dated investigation history. Older architecture notes and claims in `CLAUDE.md`, `HANDOVER.md`, `PLAN.md`, and `RESEARCH.md` can describe earlier solver versions. Consult the current source and the latest recovery report when they disagree. References to local build folders, raw logs, or external reference executables do not imply those artifacts are included in this repository.

## Build on Windows

The current application was built and launched on Windows with an NVIDIA RTX 4090 using:

| Dependency | Locally used version/setup |
| --- | --- |
| Visual Studio | 2022, MSVC C++ tools and Windows SDK |
| CMake | 4.1.1; the project declares a minimum of 3.24 |
| NVIDIA CUDA Toolkit | 13.1, including Visual Studio integration |
| Qt | 6.11.1, `msvc2022_64` package |
| Open CASCADE | 8.0, installed/built in the layout described below |
| GPU | RTX 4090, CUDA compute capability 8.9 |

The CMake configuration currently contains Windows-specific dependency discovery. Linux and macOS builds are not documented or verified.

Clone the project and run the following in PowerShell. Adjust the CUDA and Open CASCADE paths to match your machine:

```powershell
git clone https://github.com/Kulitorum/ParaCFD.git
cd ParaCFD

$env:CUDA_PATH = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$env:CUDA_PATH_V13_1 = $env:CUDA_PATH

cmake -S . -B build -G 'Visual Studio 17 2022' -A x64 `
  -DCMAKE_CUDA_ARCHITECTURES=89 `
  -DOPENCASCADE_DIR='C:/OpenCASCADE-8.0/build2' `
  -DPARACFD_ENABLE_QT=ON

cmake --build build --config Release --target paracfd-gui --parallel 8
.\build\Release\paracfd-gui.exe
```

`89` targets the RTX 4090 used for development. Set `CMAKE_CUDA_ARCHITECTURES` for your GPU; the project's configured fallback list is `75;86;89`. GPU memory requirements vary substantially with the domain, refinement, and geometry.

CMake automatically searches `C:/Qt/6.*/msvc2022_64` for Qt. For another installation, supply `-DQt6_DIR='D:/path/to/Qt/lib/cmake/Qt6'` and, if necessary, `-DWINDEPLOYQT_EXECUTABLE='D:/path/to/Qt/bin/windeployqt.exe'`. If multiple Qt versions are installed, verify that configuration and deployment select the same one.

`OPENCASCADE_DIR` currently expects headers in `inc`, import libraries in `win64/vc14/lib`, and runtime DLLs in `win64/vc14/bin`. A different installation layout requires adjusting those paths in [CMakeLists.txt](CMakeLists.txt). The post-build steps copy the Open CASCADE and Qt runtime dependencies beside the executable; retain those DLLs and plugin folders when moving it.

If CUDA detection reports an empty toolkit directory, check the two environment variables above in the shell doing the build. After a failed initial configuration, repeat the configure command with `cmake --fresh` to discard its incomplete cache.

## Run a case

Launch from the repository root so the example configuration paths resolve correctly:

```powershell
# Open the desktop application.
.\build\Release\paracfd-gui.exe

# Load and start the supplied PlanB case explicitly.
.\build\Release\paracfd-gui.exe --config configs/planb_parakite.json
```

The application can reopen the last geometry. You can also use **File → Open STEP...** or **File → Open paraglider config...**. The PlanB configuration references `Test-Data/PlanBParakite-Solid.step`; preserve its saved placement when investigating the current lift failure. STEP coordinates are interpreted in millimetres and converted to metres internally.

For work without the GUI, build the case probe:

```powershell
cmake --build build --config Release --target paraglider_case_probe --parallel 8
.\build\Release\paraglider_case_probe.exe `
  --config configs/planb_parakite.json --steps 100 --sample-every 20
```

That command is a short diagnostic run, not a lift-validation procedure. Other tools in [tools/](tools/) expose geometry checks, analytical airfoil generation, and pressure/momentum diagnostics. The scripts in `review/` are investigation tools and may require Python/NumPy, external reference programs, or adapting historical paths.

## Repository map

| Path | Contents |
| --- | --- |
| `src/core/fluid/` | AMR, CUDA flow operators, pressure and momentum updates |
| `src/core/geometry/` | STEP import, solid classification, embedded-boundary geometry |
| `src/gui/` | Desktop application, visualization, and simulation worker |
| `configs/` | Example case settings and saved placements |
| `Test-Data/` | Wing, airfoil, and smaller diagnostic geometries |
| `tools/` | Command-line probes and geometry generators |
| `review/` | Numerical investigations, scripts, saved results, and historical recovery patches |
| `tests/` | Existing diagnostic fixtures and reference material |

Generated builds, logs, simulation output, and local recovery backups are excluded from version control. Small section/panel CSV datasets in `review/` are retained as investigation evidence. Historical recovery patches are reference material; their changes may already be present in the source.

## Contributing

The most useful contribution is a reproducible explanation or repair of the lift failure. Open an issue or pull request with the geometry/configuration used, the commit and build settings, GPU, exact command or UI steps, expected physical behaviour, and observed result. Include the reasoning behind a numerical change and distinguish operator-level evidence from an end-to-end aerodynamic result.

Please preserve the supplied PlanB placement when comparing that case, and state any intentional changes to geometry, boundary conditions, momentum path, or physical parameters. Comparisons against independent methods are welcome, with their assumptions and limitations recorded.

The public handoff was built and launched successfully on 22 September 2026. The test suites were **not rerun for this handoff**. Existing reports contain historical diagnostic results; neither those results nor the successful launch establish accurate lift.

## License

ParaCFD is licensed under the **GNU General Public License, version 3 only** (`GPL-3.0-only`). See [LICENSE](LICENSE).

Third-party components retain their own licenses. The vendored [nlohmann/json](src/3rdparty/nlohmann/json.hpp) header carries its MIT license notice. Qt, Open CASCADE, and the CUDA toolkit are separate dependencies and are not relicensed by this repository.
