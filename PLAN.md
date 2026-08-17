# Paraglider CFD implementation plan

This file records actual migration status, not a claim that unfinished architecture exists.

## A. Geometry foundation — implemented

- Preserve OpenCascade STEP import and metre conversion.
- Retain source STEP face ID and optional UV data per tessellated triangle.
- Build a static triangle BVH with AABB, segment, nearest-point, and distance queries.
- Add paraglider configuration and external-aero data structures.
- Retain only mathematically reusable FP64 MAC kernels as immutable CPU/GPU validation references.

## B. Uniform brick grid — implemented foundation

- Uniform one-level block hierarchy behaves as one Cartesian grid.
- Compact brick metadata and CPU/GPU per-level coordinate hash.
- Pooled SoA host/device `Real` fields and same-level halo exchange.

The external-aero core now consumes these pooled fields directly.

## C. One-level zero-thickness EB — implemented foundation

- Exact triangle/cell clipping for a resolved planar membrane.
- Independent plus/minus fluid fragments and fabric patches.
- Split Cartesian-face aperture records and coincident blocked-face handling.
- Explicit complex-cell detection for refinement without cross-fabric merging.
- Conservative same-side small-fragment merge.

The cross-brick atlas and deterministic multi-fragment fallback are implemented. Partial Cartesian-face coverage retains separate connected apertures. Remaining: exact cross-level EB topology where fabric reaches a 2:1 interface.

## D. One-level EB projection — implemented foundation

- Matrix-free finite-volume operator with actual volume, aperture area, and distance.
- Fast implicit regular faces plus compact irregular aperture work.
- CPU-double reference/PCG and CUDA `Real` operator/CG.
- Manufactured independent-side pressure/load and projection tests.

The composite CPU reference and persistent GPU projection now cover EB divergence, RHS, pressure, regular/special-flux correction, and post-projection divergence. Disconnected active fluid components receive deterministic gauges. `ExternalAeroCore` integrates this projection into a persistent GPU timestep.

## E. Static AMR — composite pressure foundation implemented

- 2:1 hierarchy, balancing, geometry/wake refinement, hash lookup.
- Finest refinement is reserved for fabric; the configured volumetric near wake stops one level coarser.
- Halo, restriction, prolongation, and aperture-aware reflux helpers.
- Sparse cross-brick EB atlas with finest-owner selection.
- One coupled matrix-free pressure graph spanning every level.
- Compact conservative 2:1 face tiles and EB aperture connections.
- CUDA operator/composite PCG plus manufactured coupled projection tests.
- Persistent GPU divergence/RHS/correction bridge over pooled fields and compact special-flux arrays.
- Additive two-level Galerkin preconditioner with compact nonlocal aggregate connections.
- Deterministic gauges for active pressure components disconnected from the outlet.

Pending:

- aperture-aware cross-level EB topology when fabric reaches a 2:1 interface;
- recursive multigrid/V-cycle acceleration beyond the present two-level preconditioner;
- longer-run force/conservation convergence beyond the current short-case comparison.

## F. FP32 production — storage/operator foundation implemented

- New production scalar is FP32 with an optional FP64 validation build.
- Geometry preprocessing and force/reference calculations retain FP64 where valuable.
- CUDA pressure operator has measured FP32/FP64 parity and timing probes.

The bounded MacCormack/RK2 AMR advection, Smagorinsky, external-BC, and composite-projection path is GPU resident. Same-level interpolation, Smagorinsky gradients, eddy-viscosity access, and diffusion cross brick boundaries without stencil clamping; all levels build consistent forward states before correction and commit. Coarse/fine staggered-face and cell-centred samples use one-sided linear prolongation rather than a nearest-value clamp, and reverse MacCormack correction falls back to the bounded forward value when a trace changes lattice. Normal 2:1 flux tiles gather from transported fine MAC faces, then scatter projected values to the fine faces and their aperture-weighted coarse mean. Compact EB aperture velocities advance on their own same-side graph using structured carrier states whose momentum volume is the sum of adjacent composite fluid half-volumes, minmod-limited MUSCL transport on complete same-axis chains, a first-order fallback at incomplete/junction chains, molecular diffusion, and component-collocated weighted least-squares reconstruction of all nine velocity-gradient components for Smagorinsky `|S|`. Static FP32 pseudoinverses discard modes above condition number `1e4`; this prevents roundoff-amplified derivatives in nearly unsupported seam directions. Every two-sided surface patch also builds independent wall-side work records. Production reconstructs arbitrary-orientation tangential velocity, solves Spalding's smooth-wall law at the actual wall distance, applies the requested shear impulse through exact component carrier weights, and retains the equal-and-opposite CAD patch reaction. Viscous-limit, high-Re residual, dissipation, winding, and reaction-conservation gates pass in FP32/FP64; the older implicit molecular no-slip sink remains a validation twin. The combined update uses no absolute velocity cap. One global timestep is selected from a GPU maximum reduction over physical active-brick faces; ghost storage and covered bricks are excluded. Conservative same-side 0.25-volume merging plus a configurable/reportable `1e-4 h^2` numerical aperture cutoff remove compact numerical slivers. The static fabric band stores a compact six-link same-side graph and uses bounded minmod-limited linear transport wherever all required links are open, with a first-order fallback at incomplete stencils and link-restricted molecular diffusion. A separate one-level finite-volume MAC path conserves pairwise link momentum across same-level bricks and fabric blockers. Geometry-independent CPU/GPU SoA primitives extend the contract to arbitrary dual volumes and component states. Regular face volume is half of each active adjacent pressure cell; normal tiles replace the covered-side share with one quarter coarse plus one fine half-cell, and tangential links use clipped dual-face overlaps. A unified canonical topology combines both paths at refinement corners: 7,023 nodes and 7,840 transfers include 501 multi-axis states exactly once and conserve to `1.47e-14`. Its persistent GPU driver uses 6,144 retained overlap weights to derive tangential advectors from the projection's exact coarse/fine flux tiles, derives normal advectors from owned MAC states, and rebuilds 384 coarse aliases with `5.70e-8 / 4.08e-8` FP32 error (zero FP64). Alias ownership collisions fail preprocessing and 7,920 replaced directed structured links are cleared while physical boundaries remain. Production still requires compact EB states and explicit boundary-flux accounting before switching to globally conservative advection. The production smooth-wall PlanB gate remains stable for 1,500 steps / 3.831 s with a 15.57 m/s compact peak, `5.18e-6 s^-1` RMS divergence, and separated pressure/viscous/total force `[173.09,-4.73,-34.56] / [7.16,-0.024,0.033] / [180.25,-4.76,-34.53]` N. Pending: complete/refluxed composite momentum ownership and boundary fluxes plus full grid/domain/force validation.

A separate one-level external-aero finite-volume gate now applies prescribed +X inflow momentum, convective X-max export, zero Y/Z normal flux, and exact half-width boundary staggered volumes. Its uniform-freestream and localized-transverse-export manufactured cases pass in FP32 and FP64. The remaining boundary work is to attach that contract to the multilevel composite ownership graph rather than the current production semi-Lagrangian timestep.

The compact EB mass map now also has a separately callable internal first-order donor operator. It reconstructs fragment-centred vectors from aperture and regular-carrier half masses, applies equal-and-opposite `dt A u_n u_donor` transfers on embedded apertures, and scatters back to shared states. The manufactured global physical-carrier momentum error is `1.19e-7` FP32 and zero FP64. It is intentionally not a production step: an initial promotion attempt omitted the ordinary regular/compact perimeter fluxes. A deterministic ownership builder enumerates those exact open Cartesian faces with canonical MAC address, ordered pressure endpoints, area, actual composite carrier volume, and internal/perimeter flags. The persistent projection bridge uploads these addresses and evaluates the combined embedded-aperture/ordinary-MAC mass balance directly on the GPU. A finite membrane parallel to +X produces 110 records (56 true perimeter faces), closes to zero in FP32 and FP64, and reports the analytical `0.5 m^3/s` imbalance when one pooled MAC face is changed from 1 to 0.5 m/s. A separate staged GPU momentum driver now reconstructs fragment/adjacent-regular vectors from unique physical carriers, applies both aperture and ordinary-perimeter donor transfers with equal/opposite endpoint increments, scatters them once, and refreshes duplicate same-level brick-face aliases. Its perturbed 96-node/160-carrier/80-perimeter fixture conserves physical momentum to `1.80e-6` FP32 and exactly in FP64; a parallel-flow fixture preserves every compact aperture state and its 12 aliases agree exactly. The remaining promotion step is to share the perimeter flux register with the surrounding structured update so the regular endpoint's other faces close the same control-volume balance.

## G. Aerodynamics — manufactured load path implemented

- Per-triangle plus/minus pressure and Cp, delta-Cp, force, whole-wing force/moment.
- Winding-invariant pressure force.
- CL/CD/CS only with explicit positive reference area.
- Pressure drag is labelled pressure-only.

Composite patches now retain global plus/minus pressure DOFs and publish per-triangle/whole-wing pressure, smooth-wall viscous, and total loads. The GUI worker automatically transfers throttled per-triangle Cp+, Cp-, delta-Cp, and separated load summaries. Visible-side and explicit plus/minus surface rendering are implemented. Pending: establish a converged validated paraglider flow.

## H. GUI migration — implemented foundation

The executable now starts a paraglider-only Qt window. Particle/tracer fabric collision uses the triangle BVH. The viewer displays live velocity/pressure slices, arrows, tracers, static AMR/owned-EB debug boxes, and two-sided triangle pressure colours. The scalar plane samples the finest active brick at each vertex from the throttled per-level host snapshot, instead of sampling the coarse uniform tracer volume; each in-plane dimension now follows the domain extent divided by finest active `h`, subject to a defensive display-allocation cap. Cut-cell pressure samples select the corresponding fluid-fragment DOF instead of the inactive parent-cell slot, and the overlay reports the exact slice dimensions/finest spacing. Direct STEP load infers horizontal span/chord from the face-only bbox, maps the assumed negative leading-edge direction to upstream -X against the fixed +X freestream, then sizes and centres the domain around the placed wing; LE/TE polarity remains an explicit user confirmation. A shared FP32-canonical zero-origin framing helper makes GUI and probe EB clipping coordinates identical. The default close view and explicit `Fit Wing`/`Fit Domain` controls keep visual framing separate from the physically asymmetric wake domain. Successful STEP/config loads populate a persistent ten-entry recent-file menu; ordinary startup restores the last valid wing/config and placement without automatically preprocessing or running it. `Build CFD Grid + Run` performs preprocessing and immediately advances `ExternalAeroCore`; Play/Pause and Step then control it without a separate ambiguous start phase. The compact panel restores live arrow mode/density/speed/uniform-size and tracer mode/density/length/width/boring-filter controls; all editors reject wheel changes so the wheel remains dedicated to panel scrolling. Two-dimensional arrows and tracers ignore the non-depth-writing scalar plane but remain depth-tested against the STEP canopy; volume arrows retain normal occlusion. A cropped-Y diagnostic trims the open mesh at a selectable span station/requested physical width and runs the same uniform CFD path at configured finest spacing. Width snaps upward to a whole cell count (two cells / 0.125 m by default for PlanB; 1 m gives 16), the viewer selects a Y-normal 2-D view, and the overlay reports cell count and exact simulation step. The diagnostic does not alter the saved whole-wing placement/config. It reports pressure convergence/timing, separate pressure/skin-friction/total loads, regular/EB velocity, conservation, and persistent-memory diagnostics. Pending: direct CUDA/OpenGL AMR sampling (removing the throttled host copy), fragment-centred velocity reconstruction inside cut cells, and pressure/viscous-force-vector display.

## I. Obsolete-code removal — implemented

The building, centerline, roof/wall, solid voxelizer/load, channel/ground/seabed/porous, old scene/worker/window, building probe/config, old installer, and misleading legacy research paths have been deleted after the replacement release and tests passed. Git history remains the reference if an old numerical idea must be recovered.

## Definition of trustworthy version 1

- No pressure/velocity connection crosses fabric.
- Openings remain fluid-connected.
- Composite AMR flux is conservative through EB and coarse/fine interfaces.
- GPU projection reaches the declared residual tolerance with a documented norm.
- FP32 agrees with CPU/FP64 reference cases at a justified tolerance.
- Normal, parallel, and inclined dynamic plates pass; the current short equal-finest AMR inclined-plate lift differs from uniform-fine by about 0.4%.
- A missing cavity inlet admits developed bidirectional flow while its closed control has zero represented opening flux.
- Imported-wing convergence runs stop at a common physical time and emit a machine-readable summary; step count is not accepted as a proxy for elapsed flow time.
- End-to-end timestep contains no per-step CPU geometry work or bulk field transfers.
- GUI and logs clearly distinguish pressure drag from total drag and expose unresolved geometry.
