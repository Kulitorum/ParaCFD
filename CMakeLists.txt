cmake_minimum_required(VERSION 3.24)

# ScourProtection — 3D morphodynamic scour-protection simulator (C++/CUDA + Qt6).
# Milestone M0: scaffold. Three targets (libscour, scour CLI, scour-gui stub),
# GoogleTest via FetchContent, ctest wiring, JSON config loader, VTI writer.
project(ScourProtection LANGUAGES CXX CUDA)

# --- Global standards -------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA target GPUs (real SASS per arch). Default covers the common modern desktop/laptop cards so a
# shipped build runs beyond the RTX 4090 dev box: 75 = Turing (RTX 20 / GTX 16 / Quadro RTX / T4),
# 86 = Ampere (RTX 30 / A10 / A40), 89 = Ada (RTX 40, this machine). Overridable on the command line
# (e.g. -DCMAKE_CUDA_ARCHITECTURES=89 for a fast single-arch dev build). Building for 3 archs ~triples
# the device-code compile time. Older cards (Pascal 61, Volta 70) need their arch added here; for
# forward compat on unknown newer GPUs append a virtual arch (e.g. "89-virtual") to embed PTX for JIT.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
	set(CMAKE_CUDA_ARCHITECTURES 75 86 89)
endif()

# Default to Release for single-config generators (Ninja); multi-config (VS)
# selects with --config.
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

enable_testing()

# --- libscour (static core, no Qt) -----------------------------------------
add_library(libscour STATIC
	src/core/config.cpp
	src/core/vti_writer.cpp
	src/core/vti_reader.cpp
	src/core/cuda_probe.cu
	# M1 — fluid core (MAC + MacCormack + Smagorinsky + MGPCG)
	src/core/fluid/advect.cu
	src/core/fluid/turbulence.cu
	src/core/fluid/project.cu
	src/core/fluid/mgpcg.cu
	src/core/fluid/stam_fluid_core.cu
	src/core/fluid/cavity.cpp
	# M2 — open channel + obstacles (inlet/Orlanski/lid, voxel mask, cylinder gate)
	src/core/fluid/channel_ops.cu
	src/core/fluid/channel_pressure.cu
	src/core/fluid/channel_mask.cpp
	src/core/fluid/channel_core.cu
	src/core/fluid/channel_porous.cu
	src/core/fluid/cylinder.cu
	src/core/fluid/channel_empty.cpp
	# M3 — bed-shear wall function, periodic channel + body force, SEM inlet, precursor
	src/core/fluid/bedshear.cu
	src/core/fluid/channel_periodic.cu
	src/core/fluid/plane_ops.cu
	src/core/fluid/sem_inlet.cu
	src/core/fluid/channel_sem.cpp
	src/core/fluid/precursor.cpp
	# M4 — suspended sediment (c advection w/ settling, ν_t/σ_s diffusion, van Rijn pickup +
	# w_s·c_b deposition flux BC, hindered settling; V5 settling column + V6 Rouse drivers)
	src/core/sediment/suspended.cu
	src/core/sediment/sed_validation.cu
	# M5 — bed state + morphodynamics (f_pack field + states + box-filtered interface geometry,
	# Winterwerp erosion, Wong-Parker/Engelund-Fredsøe bedload donor transfers, slope-corrected
	# θ_cr, layered-column avalanching, Exner-equivalent column update, MORFAC + limiter; gate a-d)
	src/core/sediment/bedstate.cu
	src/core/sediment/avalanche.cu
	src/core/sediment/morpho_validation.cu
	# Erodible-seabed demo scenario — the FULL live morphodynamic loop (fluid→τ_b→suspended→Exner→
	# avalanche) assembled on the live flow core. seabed_morpho.cu = three bridge kernels (elevated/
	# moving-bed τ_b, suspended-field confinement, bed clamp/pin); seabed_engine.cpp = the host driver
	# that reuses the gate-verified kernels UNCHANGED. Qt-free (libscour); the GUI runs it on the worker.
	src/core/sediment/seabed_morpho.cu
	src/core/sediment/seabed_engine.cpp
	# Geometry — watertight ray-parity voxelizer (STEP TriMesh → solid obstacle mask). OCC-free
	# (needs only the mesh); the shared model placement + TriMesh headers live in core/geometry/.
	src/core/geometry/voxelize.cpp
)
set_target_properties(libscour PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_compile_features(libscour PUBLIC cxx_std_20)
# cudart_static is built against the static CRT; the rest of the build is /MD.
# Prefer the dynamic CRT (MSVCRT) everywhere to silence the benign LNK4098.
if(MSVC)
	target_link_options(libscour INTERFACE /NODEFAULTLIB:LIBCMT)
endif()
# Public: our own headers (include as "core/...") and the vendored single header.
target_include_directories(libscour PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty
)

# --- scour_geometry (STEP import via OpenCascade; isolates OCC from the core) -
# Qt-free static lib that mirrors how scour_gui_cuda splits out the GL/CUDA interop:
# it is the ONLY target that links OpenCascade, so libscour and every physics gate
# test stay OCC-free. Only scour-gui and the geometry unit test link it. STEP files
# arrive in millimetres and are scaled to metres inside step_import.cpp (PLAN §3).
set(OPENCASCADE_DIR "C:/OpenCASCADE-8.0/build2" CACHE PATH "OpenCascade 8.0 install root")
set(OCC_INCLUDE_DIR "${OPENCASCADE_DIR}/inc")
set(OCC_LIB_DIR "${OPENCASCADE_DIR}/win64/vc14/lib")
set(OCC_BIN_DIR "${OPENCASCADE_DIR}/win64/vc14/bin")

# Minimal toolkits needed to READ + TRIANGULATE a STEP shape (PLAN §3). These are the
# .lib import libs we link; being PRIVATE on a STATIC lib, CMake still propagates them
# (as $<LINK_ONLY>) to scour-gui / the geometry test so their final link resolves.
set(OCC_LINK_TOOLKITS
	TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKMesh TKXSBase TKDESTEP)

# Full transitive runtime DLL closure of the link set (computed via dumpbin /dependents:
# TKDESTEP.dll pulls in TKDE/TKXCAF/TKLCAF/TKShHealing → TKV3d/TKService → freetype, etc.).
# These + the 3rdparty DLLs must sit next to any exe that loads scour_geometry so it stays
# double-clickable with NO OpenCascade on PATH (see CLAUDE.md).
set(OCC_RUNTIME_TOOLKITS
	TKBO TKBRep TKCAF TKCDF TKDE TKDESTEP TKernel TKG2d TKG3d TKGeomAlgo TKGeomBase
	TKHLR TKLCAF TKMath TKMesh TKPrim TKService TKShHealing TKTopAlgo TKV3d TKVCAF TKXCAF TKXSBase)
set(OCC_3RDPARTY_DLLS freetype brotlicommon brotlidec bz2 libpng16 zlib1)

# Copy the OpenCascade runtime DLLs next to a target's exe (POST_BUILD, copy_if_different so
# it is a no-op once deployed). Called for both scour-gui and the geometry test exe.
function(scour_deploy_occ_dlls target)
	if(NOT WIN32)
		return()
	endif()
	set(_paths)
	foreach(_tk ${OCC_RUNTIME_TOOLKITS} ${OCC_3RDPARTY_DLLS})
		set(_p "${OCC_BIN_DIR}/${_tk}.dll")
		if(EXISTS "${_p}")
			list(APPEND _paths "${_p}")
		else()
			message(WARNING "scour_deploy_occ_dlls(${target}): OCC DLL not found: ${_p}")
		endif()
	endforeach()
	add_custom_command(TARGET ${target} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_paths} "$<TARGET_FILE_DIR:${target}>"
		COMMENT "Deploying OpenCascade runtime DLLs next to ${target} (double-clickable, no OCC on PATH)"
		VERBATIM)
endfunction()

add_library(scour_geometry STATIC src/core/geometry/step_import.cpp)
set_target_properties(scour_geometry PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(scour_geometry PUBLIC cxx_std_20)
# OCC headers are PRIVATE — consumers include only the OCC-free step_import.h.
target_include_directories(scour_geometry
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
	PRIVATE ${OCC_INCLUDE_DIR})
# OCCT 8.0 moved some typedefs into a Deprecated/NCollectionAliases dir; include it +
# silence the deprecation notices (matches the OCC 8.0 usage in cobod-slicer).
set(OCC_DEPRECATED_DIR "${OPENCASCADE_DIR}/../OCCT/src/Deprecated/NCollectionAliases")
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(scour_geometry PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(scour_geometry PRIVATE OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
foreach(_tk ${OCC_LINK_TOOLKITS})
	if(NOT EXISTS "${OCC_LIB_DIR}/${_tk}.lib")
		message(FATAL_ERROR "scour_geometry: OpenCascade import lib not found: ${OCC_LIB_DIR}/${_tk}.lib")
	endif()
	target_link_libraries(scour_geometry PRIVATE "${OCC_LIB_DIR}/${_tk}.lib")
endforeach()

# --- scour CLI --------------------------------------------------------------
add_executable(scour src/cli/main.cpp)
target_link_libraries(scour PRIVATE libscour)

# --- scour_gui_cuda (GUI CUDA half: slice colour kernel + CUDA-GL interop) ---
# Kept Qt-free so nvcc never sees Qt's host flags, and so the fast unit suite can
# link the slice kernel + its CPU reference without pulling in Qt/GL.
add_library(scour_gui_cuda STATIC
	src/gui/slice_field.cu   # pure CUDA sampler+colourmap (+ CPU reference)
	src/gui/slice_gl.cu      # cudaGraphicsGLRegisterBuffer interop (the only GL-header TU)
)
set_target_properties(scour_gui_cuda PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(scour_gui_cuda PUBLIC libscour)

# --- scour-gui (Qt6 + GL 4.3 slice viewer; G1) ------------------------------
# Default ON from G1 onward (Qt 6.11.1 present on this machine). Turn OFF to build
# a Qt-less stub (main_stub.cpp) on hosts without Qt.
option(SCOUR_ENABLE_QT "Build scour-gui against Qt6 (G1 slice viewer)" ON)
add_executable(scour-gui)
if(SCOUR_ENABLE_QT)
	file(GLOB _qt_candidates "C:/Qt/6.*/msvc2022_64")
	list(SORT _qt_candidates)
	list(REVERSE _qt_candidates) # newest first
	foreach(_qt_dir ${_qt_candidates})
		if(EXISTS "${_qt_dir}/lib/cmake/Qt6/Qt6Config.cmake")
			list(PREPEND CMAKE_PREFIX_PATH "${_qt_dir}")
			set(_qt_found "${_qt_dir}")
			break()
		endif()
	endforeach()
	find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets)
	message(STATUS "scour-gui: Qt6 from ${_qt_found}")
	target_sources(scour-gui PRIVATE
		src/gui/main.cpp
		src/gui/main_window.cpp
		src/gui/slice_viewer.cpp
		src/gui/sim_worker.cpp
		src/gui/sim_setup.cpp
		src/gui/scene_io.cpp
		src/gui/flow_particles.cpp
		src/gui/flow_tracers.cpp
		src/gui/video_recorder.cpp
		src/gui/video_settings_dialog.cpp
	)
	set_target_properties(scour-gui PROPERTIES AUTOMOC ON)
	target_link_libraries(scour-gui PRIVATE
		libscour scour_gui_cuda scour_geometry
		Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
	)
	target_compile_definitions(scour-gui PRIVATE SCOUR_HAVE_QT=1)

	# scour_geometry links OpenCascade DLLs -> deploy them next to scour-gui.exe so the
	# exe stays double-clickable (windeployqt only handles Qt, not OCC).
	scour_deploy_occ_dlls(scour-gui)

	# --- Deploy the Qt runtime next to the exe so scour-gui.exe is double-clickable
	# (no Qt on PATH needed) and stays so after a clean rebuild. windeployqt copies the
	# Qt6 DLLs + the platforms/qwindows.dll plugin into the output dir. Runs after each
	# link; it is a no-op when the deployed files are already up to date.
	find_program(WINDEPLOYQT_EXECUTABLE windeployqt HINTS "${_qt_found}/bin")
	if(WINDEPLOYQT_EXECUTABLE)
		add_custom_command(TARGET scour-gui POST_BUILD
			COMMAND "${WINDEPLOYQT_EXECUTABLE}"
				"$<IF:$<CONFIG:Debug>,--debug,--release>"
				--no-translations --no-system-d3d-compiler "$<TARGET_FILE:scour-gui>"
			COMMENT "windeployqt: deploying Qt runtime next to scour-gui.exe (double-clickable)"
			VERBATIM)
		# windeployqt only ships the windows platform plugin; also deploy qoffscreen so
		# the headless --offscreen smoke path (STEP-load CI) works with no Qt on PATH.
		set(_qoffscreen "${_qt_found}/plugins/platforms/$<IF:$<CONFIG:Debug>,qoffscreend,qoffscreen>.dll")
		add_custom_command(TARGET scour-gui POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:scour-gui>/platforms"
			COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_qoffscreen}" "$<TARGET_FILE_DIR:scour-gui>/platforms"
			COMMENT "Deploying Qt offscreen platform plugin next to scour-gui.exe"
			VERBATIM)
	else()
		message(WARNING "scour-gui: windeployqt not found under ${_qt_found}/bin; "
			"the exe will need the Qt bin dir on PATH to run (see CLAUDE.md).")
	endif()

	# --- G1 automated gate: launch, run the bundled config offscreen-capable for
	# ~10 s, exit 0 on scripted close (manual fps + GL-thread assert per PLAN §4).
	add_test(NAME gate_G1
		COMMAND scour-gui --config "${CMAKE_CURRENT_SOURCE_DIR}/configs/g1_viewer.json" --autoclose-ms 10000
	)
	set_tests_properties(gate_G1 PROPERTIES LABELS gate_G1 TIMEOUT 120)
else()
	target_sources(scour-gui PRIVATE src/gui/main_stub.cpp)
	target_link_libraries(scour-gui PRIVATE libscour)
endif()

# --- GoogleTest via FetchContent -------------------------------------------
include(FetchContent)
FetchContent_Declare(googletest
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG        v1.15.2
	GIT_SHALLOW    TRUE
)
# Match MSVC dynamic CRT (/MD) used by the rest of the build.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(scour_tests
	tests/test_config.cpp
	tests/test_vti.cpp
	tests/test_cuda.cpp
	# M1 — GPU-vs-CPU kernel parity + MGPCG solver tests (CUDA sources)
	tests/test_fluid_kernels.cu
	tests/test_mgpcg.cu
	# M2 — masked open-channel kernel parity + ChannelMgpcg convergence
	tests/test_channel_kernels.cu
	# M3 — bed-shear wall function, periodic-channel closure, SEM, precursor
	tests/test_bedshear.cu
	tests/test_m3.cu
	# M4 — suspended-sediment kernel parity + sed_physics value checks
	tests/test_suspended.cu
	# M5 — bed state + morphodynamics kernel parity + sed_physics M5 value checks + avalanche
	tests/test_bedstate.cu
	# Seabed demo — GPU-vs-CPU parity for the three seabed bridge kernels (τ_b, confine, clamp/pin)
	tests/test_seabed.cu
	# Geometry — ray-parity voxelizer volume + inside/outside gates (OCC-free, procedural cube)
	tests/test_voxelize.cpp
	# G1 — slice colour kernel GPU-vs-CPU parity (gui/slice_field.cu)
	tests/test_slice_field.cu
	# M7 (M-DEP) — trapping-efficiency engineering cross-checks (van Rijn e_s, Camp-Hazen, screen k)
	tests/test_deposition.cpp
	# M7 — sub-grid porous momentum-sink GPU-vs-CPU parity + drag limits (channel_porous.cu)
	tests/test_porous.cu
	# M8 — campaign machinery: timescale/MORFAC/equilibrium-fit/ranking maths (run_protocol.h)
	tests/test_campaign.cpp
)
set_target_properties(scour_tests PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_tests PRIVATE libscour scour_gui_cuda GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(scour_tests PROPERTIES LABELS unit DISCOVERY_TIMEOUT 120)

# --- scour_geometry_tests (STEP import unit test; links OCC, kept separate) --
# A small exe so the physics scour_tests target stays OCC-free. Loads the bundled
# tests/inputs/cube.step and checks the mm->m scaling + triangle count. DISCOVERY_MODE
# PRE_TEST defers test enumeration to ctest time (after POST_BUILD deploys the OCC DLLs
# next to the exe), so discovery does not race the DLL copy.
add_executable(scour_geometry_tests tests/test_step_import.cpp)
target_link_libraries(scour_geometry_tests PRIVATE scour_geometry GTest::gtest GTest::gtest_main)
target_compile_definitions(scour_geometry_tests PRIVATE
	STEP_CUBE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/inputs/cube.step")
scour_deploy_occ_dlls(scour_geometry_tests)
gtest_discover_tests(scour_geometry_tests PROPERTIES LABELS unit DISCOVERY_TIMEOUT 120 DISCOVERY_MODE PRE_TEST)

# --- scour_physics (Jolt rigid-body drop/settle prep phase; isolates Jolt from the core) ----
# The ONLY target that links Jolt Physics — exactly as scour_geometry is the only one that links
# OpenCascade — so libscour and every physics gate stay Jolt-free. Settling runs once on the CPU
# and emits ModelPlacements consumed by the existing voxelizer (PLAN G3). Uses the prebuilt
# Release /MD install (find_package CONFIG imports Jolt::Jolt with the AVX2 flag + JPH_* defines
# auto-applied, so this consumer stays define-matched). Optional so hosts without Jolt still build.
option(SCOUR_ENABLE_JOLT "Build the PLAN G3 drop/settle preparation phase (Jolt Physics)" ON)
if(SCOUR_ENABLE_JOLT)
	set(JOLT_INSTALL_DIR "C:/CODE/JoltInstall" CACHE PATH "Prebuilt Jolt Physics install root (find_package CONFIG)")
	find_package(Jolt CONFIG QUIET PATHS "${JOLT_INSTALL_DIR}/lib/cmake/Jolt")
	if(Jolt_FOUND)
		message(STATUS "scour_physics: Jolt from ${JOLT_INSTALL_DIR} (target Jolt::Jolt)")
		add_library(scour_physics STATIC src/core/physics/settle.cpp)
		set_target_properties(scour_physics PROPERTIES POSITION_INDEPENDENT_CODE ON)
		target_compile_features(scour_physics PUBLIC cxx_std_20)
		# Jolt is PRIVATE — consumers include only the Jolt-free settle.h. Being PRIVATE on a STATIC
		# lib, CMake still propagates Jolt.lib ($<LINK_ONLY>) to the final link (same as OCC above).
		target_include_directories(scour_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
		target_link_libraries(scour_physics PRIVATE Jolt::Jolt)

		# G3.1 gate (fast — runs in the unit suite): drop a scatter of blocks, settle, and assert
		# they sleep, rest on the bed, don't interpenetrate, and reproduce bit-for-bit for a seed.
		add_executable(scour_physics_tests tests/test_settle.cpp)
		target_link_libraries(scour_physics_tests PRIVATE scour_physics GTest::gtest GTest::gtest_main)
		gtest_discover_tests(scour_physics_tests PROPERTIES LABELS unit DISCOVERY_TIMEOUT 120)

		# Wire the drop/settle preparation phase into the GUI (PLAN G3.2). SCOUR_HAVE_JOLT guards the
		# drop code so scour-gui still compiles on a host without Jolt (the feature is then absent).
		if(TARGET scour-gui)
			target_link_libraries(scour-gui PRIVATE scour_physics)
			target_compile_definitions(scour-gui PRIVATE SCOUR_HAVE_JOLT=1)
			message(STATUS "scour-gui: drop/settle preparation phase enabled (Jolt)")
		endif()
	else()
		message(WARNING "scour_physics: Jolt not found at ${JOLT_INSTALL_DIR}; the PLAN G3 "
			"drop/settle prep phase is SKIPPED. Set -DJOLT_INSTALL_DIR=<install> to enable.")
	endif()
endif()

# --- M1 gate: gate_M1 -------------------------------------------------------
# Lid-driven cavity vs Ghia (1982). Long-running (minutes on the RTX 4090):
# kept out of the fast unit suite as its own executable + ctest.
add_executable(scour_m1_gate tests/gate_m1_main.cpp)
set_target_properties(scour_m1_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m1_gate PRIVATE libscour)
add_test(NAME gate_M1
	COMMAND scour_m1_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/v1_cavity_re100.json"
)
set_tests_properties(gate_M1 PROPERTIES LABELS gate_M1 TIMEOUT 1800)

# --- M2 gate: gate_M2 -------------------------------------------------------
# Cylinder Cd/St (Re=100 & 200) + open-channel BCs vs research/11 §4. Long-running
# (many minutes on the RTX 4090): its own executable + ctest.
add_executable(scour_m2_gate tests/gate_m2_main.cpp)
set_target_properties(scour_m2_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m2_gate PRIVATE libscour)
add_test(NAME gate_M2
	COMMAND scour_m2_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/v3_cylinder.json"
)
set_tests_properties(gate_M2 PROPERTIES LABELS gate_M2 TIMEOUT 5400)

# --- M3 gate: gate_M3 -------------------------------------------------------
# Flat-bed periodic channel log-law (bed-shear wall function) + Jarrin SEM inlet vs
# RESEARCH §4 anchor / §8. Its own executable + ctest.
add_executable(scour_m3_gate tests/gate_m3_main.cpp)
set_target_properties(scour_m3_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m3_gate PRIVATE libscour)
add_test(NAME gate_M3
	COMMAND scour_m3_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/v2_channel_loglaw.json"
)
set_tests_properties(gate_M3 PROPERTIES LABELS gate_M3 TIMEOUT 5400)

# --- M4 gate: gate_M4 -------------------------------------------------------
# Suspended sediment: V5 still-water settling column (L2 + mass conservation) and V6
# equilibrium Rouse profile (fitted exponent vs R=w_s/(β·κ·u*)). RESEARCH §6.1-6.2. The
# tiny validation grids run in seconds, but it is kept out of the fast unit suite as a gate.
add_executable(scour_m4_gate tests/gate_m4_main.cpp)
set_target_properties(scour_m4_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m4_gate PRIVATE libscour)
add_test(NAME gate_M4
	COMMAND scour_m4_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m4_sediment.json"
)
set_tests_properties(gate_M4 PROPERTIES LABELS gate_M4 TIMEOUT 900)

# --- M5 gate: gate_M5 -------------------------------------------------------
# Bed state + morphodynamics: (a-i) still-water deposition column, (a-ii) prescribed bedload vs
# 1D Exner, (b) sand-pile avalanche, (c) three-reservoir mass conservation over 1e4 steps,
# (d) Soulsby-Whitehouse threshold onset. RESEARCH §5, §7. The 1e4-step (c) run on the RTX 4090
# takes a few minutes; its own executable + ctest, kept out of the fast unit suite.
add_executable(scour_m5_gate tests/gate_m5_main.cpp)
set_target_properties(scour_m5_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m5_gate PRIVATE libscour)
add_test(NAME gate_M5
	COMMAND scour_m5_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m5_morpho.json"
)
set_tests_properties(gate_M5 PROPERTIES LABELS gate_M5 TIMEOUT 1800)

# --- M6 gate: gate_M6 -------------------------------------------------------
# Pile-scour calibration (Roulund et al. 2005 / RESEARCH V7): the live morphodynamic loop around a
# voxelized vertical pile in an erodible sand bed. Frozen-bed u* pre-check → MORFAC-ramped morphology
# → S(t) dual-fit → S/D vs the V7 band. SUPERVISED: the gate REPORTS S/D + sweeps the calibration
# knobs (α, Cs, Winterwerp erosion_coeff, HSV) but never freezes a constant (CLAUDE.md M6). The full
# 2.9M-cell M=1 run takes many minutes on the RTX 4090 → its own executable + ctest (background+poll).
add_executable(scour_m6_gate tests/gate_m6_main.cpp)
set_target_properties(scour_m6_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m6_gate PRIVATE libscour)
add_test(NAME gate_M6
	COMMAND scour_m6_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m6_roulund.json"
)
set_tests_properties(gate_M6 PROPERTIES LABELS gate_M6 TIMEOUT 7200)

# --- M7 gate: gate_M7 (M-DEP / Du 2025 / V8) --------------------------------
# Deposition/protection gate: bare pile vs C-AR / H-AR perforated porous-unit arrays (sub-grid
# thin-screen model). Comparative/qualitative: porous units must reduce scour + trap sand; absolute
# S/D calibration is supervised (gate_M6 report). Three ~2.9M-cell cases → many minutes; background+poll.
add_executable(scour_m7_gate tests/gate_m7_main.cpp)
set_target_properties(scour_m7_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m7_gate PRIVATE libscour)
add_test(NAME gate_M7
	COMMAND scour_m7_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m7_du.json"
)
set_tests_properties(gate_M7 PROPERTIES LABELS gate_M7 TIMEOUT 7200)

# --- M8 gate: gate_M8 (shape ranking) ---------------------------------------
# Runs procedural reference shapes (plate/ring/cup) on the live pipeline, fits trapped-volume
# equilibria with CIs, and ranks. Comparative gate: cup out-traps plate + a distinguishable pair.
# Compact ranking patch → a few minutes; background+poll.
add_executable(scour_m8_gate tests/gate_m8_main.cpp)
set_target_properties(scour_m8_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m8_gate PRIVATE libscour)
add_test(NAME gate_M8
	COMMAND scour_m8_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m8_ranking.json"
)
set_tests_properties(gate_M8 PROPERTIES LABELS gate_M8 TIMEOUT 3600)

# --- M9 gate: gate_M9 (tidal reversal + ranking campaign) --------------------
# Face-swap tidal reversal (RESEARCH §8, cosine slack ramps, MORFAC ≤ 5): checks total-sand
# conservation < 0.1 % (CLOSED boundary ⇒ machine precision), no monotonic KE growth across the 4
# slack transitions, and the MORFAC-limiter clip-fraction < 1 %; then ranks ≥ 3 candidate shapes by
# matched-time trapped volume under the identical reversing tide. Compact patch → minutes; bg+poll.
add_executable(scour_m9_gate tests/gate_m9_main.cpp)
set_target_properties(scour_m9_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_m9_gate PRIVATE libscour)
add_test(NAME gate_M9
	COMMAND scour_m9_gate "${CMAKE_CURRENT_SOURCE_DIR}/configs/m9_tidal.json"
)
set_tests_properties(gate_M9 PROPERTIES LABELS gate_M9 TIMEOUT 5400)

# --- Voxelizer flow-integration gate: voxel_flow ----------------------------
# Headless end-to-end check that a voxelized STEP obstacle actually diverts the live M2 flow:
# a procedural box is voxelized into a small open channel, stepped a few hundred steps, and the
# settled field is asserted (no penetration into solid cells; peak fluid speed > inlet U). Links
# only libscour (OCC-free, Qt-free); runs in seconds on the RTX 4090 as its own labelled gate.
add_executable(scour_voxel_flow_gate tests/gate_voxel_flow_main.cpp)
set_target_properties(scour_voxel_flow_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_voxel_flow_gate PRIVATE libscour)
add_test(NAME voxel_flow COMMAND scour_voxel_flow_gate)
set_tests_properties(voxel_flow PROPERTIES LABELS voxel_flow TIMEOUT 600)

# --- Seabed morphodynamic-loop verification gate: seabed_flat -----------------
# Headless (no Qt/GL) end-to-end check that the live fluid→τ_b→suspended→Exner→avalanche loop is
# STABLE and MASS-CONSERVING: a flat 1 m sand bed under a current (no structure) and a variant with
# a procedural box structure are each stepped a few hundred fluid+morpho steps. Asserts total sand
# (bed + suspended) is conserved, the flow stays bounded (no blow-up), and — with the box — the bed
# LOWERS near the obstacle (scour signature). Links only libscour; runs in ~1 min on the RTX 4090.
add_executable(scour_seabed_gate tests/gate_seabed_main.cpp)
set_target_properties(scour_seabed_gate PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
target_link_libraries(scour_seabed_gate PRIVATE libscour)
add_test(NAME seabed_flat COMMAND scour_seabed_gate)
set_tests_properties(seabed_flat PROPERTIES LABELS seabed_flat TIMEOUT 900)

# --- Gate: gate_M0 ----------------------------------------------------------
# Producer/consumer round-trip: the CLI writes a synthetic-field .vti, then the
# C++ re-parser and (if present) the Python validator confirm header + payload.
set(M0_VTI "${CMAKE_BINARY_DIR}/m0_synthetic.vti")

add_test(NAME gate_M0_emit COMMAND scour --emit-test-vti "${M0_VTI}")
set_tests_properties(gate_M0_emit PROPERTIES
	LABELS gate_M0
	FIXTURES_SETUP M0_VTI_FIXTURE
)

# C++ validator (always available): re-parse and compare against the analytic field.
add_test(NAME gate_M0_validate_cpp COMMAND scour_tests --gtest_filter=VtiRoundTrip.*)
set_tests_properties(gate_M0_validate_cpp PROPERTIES LABELS gate_M0)

# Python validator (independent re-implementation), only if an interpreter exists.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
	add_test(NAME gate_M0_validate_py
		COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_vti.py" "${M0_VTI}" --check-synthetic
	)
	set_tests_properties(gate_M0_validate_py PROPERTIES
		LABELS gate_M0
		FIXTURES_REQUIRED M0_VTI_FIXTURE
	)
	message(STATUS "gate_M0: Python validator enabled (${Python3_EXECUTABLE})")
else()
	message(STATUS "gate_M0: Python not found — C++ validator only")
endif()
