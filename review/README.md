Latest recovery results: [Final attempt — pressure traction, small-cell increments and independent checks](FINAL_ATTEMPT.md). The real-flight lift remains unresolved.

Current implementation and validation results: [REPAIR_STATUS.md](REPAIR_STATUS.md). The audit below describes the pre-repair state.

# Reproducing the September 2026 review

The review adds diagnostics, not a solver fix. [SIMULATION_REVIEW.md](../SIMULATION_REVIEW.md) records findings and measured results; [COMPLETION_PLAN.md](../COMPLETION_PLAN.md) gives acceptance gates. The source snapshot included pre-existing uncommitted work. [source-fingerprints.json](source-fingerprints.json) identifies the principal tested inputs and numerical files.

Build the current libraries/probes from the repository root on the configured Windows/CUDA/OCCT machine:

```powershell
$env:CUDA_PATH='C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$env:CudaToolkitDir=$env:CUDA_PATH
cmake --build build-paraglider-ui --config Release --target paraglider_case_probe naca_validation_probe paraglider_gpu_probe parity_probe aero_convergence_probe visualization_probe -- /m /p:TrackFileAccess=false
ctest --test-dir build-paraglider-ui -C Release --output-on-failure
```

Build and run the standalone review diagnostic against those libraries:

```powershell
.\review\build-audit.ps1
.\build-paraglider-ui\Release\physics_audit_probe.exe
.\build-paraglider-ui\Release\physics_audit_probe.exe configs/naca_validation.json --naca-slab
```

The script uses the local VS2022/CUDA13.1/OCCT8 installation and puts executable/object output in the existing ignored build directory. It does not edit CMake or production sources. The audit exits nonzero on execution failure, but **reports numerical defects rather than asserting they pass**. Its zero-time trials reset the initial/evolved fields between calls, use zero molecular/SGS viscosity, and do not apply the wall sink. The manufactured pressure-reaction fixture has zero pressure on the external boundaries. The optional slab fixture uses zero incidence and performs topology diagnostics only; it is distinct from the +4° time-limited benchmark below.

Commands executed during the review (historical record):

```powershell
.\build-paraglider-ui\Release\paraglider_case_probe.exe --config configs/planb_parakite.json --physical-time 1.6 --sample-every 100 --staggered-momentum
.\build-paraglider-ui\Release\paraglider_case_probe.exe --config configs/planb_parakite.json --physical-time 1.6 --sample-every 200 --aoa 12 --staggered-momentum
.\build-paraglider-ui\Release\paraglider_case_probe.exe --config configs/planb_parakite.json --steps 100 --sample-every 20 --aoa 12 --conservative-cell-momentum --momentum-balance
.\build-paraglider-ui\Release\naca_validation_probe.exe --profile 2412 --single-angle 4 --max-wall-seconds 180 --progress-steps 500
```

The NACA run deliberately reported a controlled timeout, not a completed coefficient gate. PlanB endpoint forces are instantaneous, and the collocated and MAC runs cover different physical durations. GPU run performance depends on the workstation load.

The user subsequently confirmed that PlanB's supplied AoA is the neutral line trim and must remain unchanged, with approximately 10% sink expected. The two `--aoa 12` commands above were earlier temporary diagnostics; they never wrote the configuration and should not be repeated as part of the recovery acceptance case. Use the unchanged-trim command and the default standalone operator audit for follow-up work. Angles on separate NACA benchmark geometries are independent of PlanB's trim.

Local raw evidence (ignored by the existing `*.log` rule): `ctest.log`, `physics-audit.log`, `naca-solid-mask-audit.log`, `planb-mac-0deg.log`, `planb-mac-12deg.log`, `planb-cv-12deg.log`, and `naca2412-mac-4deg.log`. Important measurements are preserved in the review document so they remain available without these logs. Build output is in `../review-build.log` and `audit-build.log`.
# Latest flight-recovery follow-up

See [FLIGHT_RECOVERY.md](FLIGHT_RECOVERY.md) for the latest numerical fixes, accepted and rejected runs, and the unresolved flight-condition metadata. The supplied PlanB trim remains unchanged, and a carried-weight flight result is not yet validated.
