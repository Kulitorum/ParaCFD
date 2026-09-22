Latest recovery results: [Final attempt — pressure traction, small-cell increments and independent checks](FINAL_ATTEMPT.md). The real-flight lift remains unresolved.

# Flight recovery follow-up — 6 September 2026

The flight result remains **unresolved**. This follow-up fixes additional numerical consistency defects; it does not claim that the simulated wing supports a person. The supplied PlanB CAD, placement matrix, neutral trim and case JSON are unchanged.

## Additional defects and repairs

The cell-to-face velocity reconstruction previously interpolated between cell centres using only the face-normal coordinate. A cut aperture or AMR face centroid generally lies off that line. The transfer therefore changed even a known affine velocity field. The new production test failed on both CPU and GPU with a maximum error of **0.0160576 m/s** on the inclined-box fixture.

The first full-gradient correction passed the small fixture but failed on the actual PlanB geometry: its one-ring velocity gradient stencil lost a direction on some thin fragments, leaving an affine error of **0.0009266 m/s**. An evolved PlanB run also developed unbounded velocity near small surface fragments. Its positive force readings are invalid. That run was stopped and is retained in `flight-planb-affine-h125.log`.

The final reconstruction uses the pressure topology's checked, same-fluid gradient support, which can extend beyond one ring. It adds the correction from the line between cell centres to the actual face centroid. A minmod limiter retains the correction supported by both adjacent gradients. Identical affine gradients reproduce the exact face value; disagreement limits the extrapolation. CPU and GPU use the same method. The obsolete, unused `restore_projected_fluxes` API was removed because it assumed the previous interpolation formula.

The new tests also exercise cell momentum with zero mass flux, zero pressure and vanishing diffusion/wall timesteps. These supplement the previous tiny-timestep test, which only exercised MAC aperture transport even when initialized from a cell-momentum run. They do not establish the accuracy of every small-cell stabilization configuration: the tested full PlanB coarse case has no momentum-redistribution groups. The remaining full-state averaging in configurations with such groups still needs a dedicated affine-preservation repair/validation.

## Evidence

The final Release build, including the GUI, succeeds. All eight CTests pass (5.51 s); the separate full-PlanB operator check also passes. The final build only removes an unused diagnostic API after the numerical runs below; it does not alter their executed solver path. Interactive GUI lifecycle testing and an independent 3D reference comparison remain outstanding.

| Check | Result |
|---|---:|
| Full PlanB affine face-velocity error, GPU | 6.4373e-7 m/s |
| Full PlanB affine face-velocity error, CPU | 5.2082e-7 m/s |
| Full PlanB pressure/body impulse mismatch | 2.9617e-8 relative |
| Full PlanB constant-pressure derivative | 0 |
| Full PlanB linear-pressure derivative, worst axis | 2.9802e-7 |
| Cell no-op/tiny-timestep affine changes on coarse PlanB | 0 at stored precision |
| NACA2412 +4°, inviscid, 128 cells/chord | Existing gates PASS |
| NACA mean CL / recorded reference | 0.71030 / 0.7414 |
| NACA force-window drift / circulation CL | 0.280% / 0.64433 |
| Final unchanged PlanB regression | 1.20032 s, 1,245 steps, bounded velocity |
| Final PlanB force in current horizontal flow | [43.1603, 0.1105, -37.9503] N |

The accepted airfoil result remains one inviscid benchmark, with approximately 4.2% lift error. It is not independent validation of the complete 3D, viscous paraglider. The PlanB result is an endpoint regression value, not a flight-load average with a mesh/time uncertainty estimate.

A finer **pre-follow-up-reconstruction** PlanB run at h=0.0625 m completed 1.20014 s / 3,460 steps with force [29.4374, 0.0317, -20.9209] N. The large difference from the earlier h=0.125 m result demonstrates unresolved mesh sensitivity. This is not a convergence result for the final reconstruction: the finer case must be rerun with the final method before forming such a comparison.

The exact, closed midspan triangle intersection is now exported as segments as well as the older vertical upper/lower envelope. Its 800 segments form one closed contour, with every endpoint incident to two segments. Independent XFOIL attempts on the exact contour, with multiple panel counts and direct panel copying, still produced large, nonconvergent spurious pressure drag. These results are rejected. They cannot be used to decide the sign or magnitude of PlanB lift. The new section export does not change the flow geometry.

## Flight condition that must be resolved

The current product fixes the wing in **horizontal +X airflow**, with zero normal flow at the upper/lower and lateral domain boundaries. It does not simulate a descending aircraft. Neutral rigging fixes the wing/line configuration; the aerodynamic incidence also depends on the direction of relative airflow. [NASA's definition](https://www.grc.nasa.gov/www/k-12/Aero2000/studweb/2-3-2i2.html) relates angle of attack to relative wind, rather than to the horizon.

For 10% sink in still air, a horizon-frame wing moving horizontally at 10 m/s and descending at 1 m/s sees a 1 m/s upward component of relative airflow. Its inclination is atan(0.1), approximately 5.71°. This is a kinematic consequence, not a lift correction or a proposed change to neutral rigging. If the exported geometry is already aligned to the incoming air, adding that inclination again would change the user's supplied AoA. **The horizon-versus-incoming-air export convention has been asked and remains unanswered. No convention or airflow offset has been assumed.**

Once that convention is established, the implementation path is:

1. Preserve the neutral CAD placement separately from flight velocity and gravity. Represent the actual relative-airflow vector, with consistent open inflow/outflow conditions on the relevant domain faces. Merely adding a vertical inlet velocity while retaining impermeable upper/lower boundaries is incorrect.
2. Obtain the measured neutral-trim airspeed, all-up mass and flight conditions. Match that case, including line/pilot/harness drag where needed. Check force balance for the complete aircraft; do not prescribe lift or tune an AoA to produce it.
3. Run the unchanged-trim case through three meshes, smaller timesteps, longer averaging windows and domain-size checks, plus the planned independent finite-wing comparison. Accept a carried-weight/glide prediction only when those checks and the defined flight measurements agree within a reported uncertainty.

## Artifacts and reproduction

`flight-recovery.patch` records only this follow-up's code changes relative to `build-paraglider-ui/flight-repair-baseline/`. `flight-fingerprints.json` records the final source and binary hashes. Earlier repairs and limitations remain documented in [REPAIR_STATUS.md](REPAIR_STATUS.md) and [../COMPLETION_PLAN.md](../COMPLETION_PLAN.md).

Raw logs are retained under `review/flight-*.log`. Runs named `limited` were stopped after the actual-geometry affine test identified the separate stencil defect; they are not accepted results. The accepted final numerical runs are `flight-planb-checked-contract.log`, `flight-planb-checked-h125.log` and `flight-naca2412-checked.log`.

```powershell
$env:CUDA_PATH = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1'
$env:CudaToolkitDir = $env:CUDA_PATH
cmake --build build-paraglider-ui --config Release -- /m
ctest --test-dir build-paraglider-ui -C Release --output-on-failure
$env:PARACFD_EB_THREADS = '16'
.\build-paraglider-ui\Release\physics_contract_probe.exe --planb
.\build-paraglider-ui\Release\paraglider_case_probe.exe --config configs/planb_parakite.json --physical-time 1.2 --sample-every 100 --projection-tolerance 1e-5
.\build-paraglider-ui\Release\naca_validation_probe.exe --profile 2412 --single-angle 4 --max-levels 3 --span-layers 4 --physical-time 1.2 --conservative-cell-momentum --max-wall-seconds 600 --progress-steps 1000
```
