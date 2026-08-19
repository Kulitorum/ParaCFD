cmake_minimum_required(VERSION 3.24)

# ParaCFD — GPU paraglider aerodynamics from zero-thickness STEP fabric
# (C++20/CUDA + Qt6/OpenGL + OpenCascade).
project(ParaCFD LANGUAGES CXX CUDA)

option(PARACFD_VALIDATION_FP64 "Use FP64 fields in the new paraglider solver path" OFF)

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

# CUDA toolkit headers (cuda_runtime.h) for the host .cpp TUs that include them.
# nvcc adds these implicitly for .cu sources, but MSVC compiling our .cpp needs the
# path explicitly; exported PUBLIC from libparacfd so all consumers pick it up.
find_package(CUDAToolkit REQUIRED)

# --- libparacfd (static core, no Qt) -----------------------------------------
add_library(libparacfd STATIC
	src/core/paraglider_config.cpp
	src/core/aero_loads.cpp
	src/core/cuda_probe.cu
	# Reusable FP64 reference kernels retained for CPU/GPU parity validation.
	src/core/fluid/advect.cu
	src/core/fluid/turbulence.cu
	src/core/fluid/project.cu
	src/core/fluid/mgpcg.cu
	# Production paraglider path: static 2:1 bricks and pooled FP32 SoA fields.
	src/core/fluid/amr_grid.cpp
	src/core/fluid/amr_fields.cu
	src/core/fluid/amr_advection.cu
	src/core/fluid/amr_exchange.cpp
	src/core/fluid/amr_eb.cpp
	src/core/fluid/amr_eb.cu
	src/core/fluid/amr_pressure.cpp
	src/core/fluid/amr_pressure.cu
	src/core/fluid/external_aero_core.cpp
	src/core/fluid/eb_pressure.cpp
	src/core/fluid/eb_pressure.cu
	# Paraglider geometry acceleration and zero-thickness embedded boundaries.
	src/core/geometry/triangle_bvh.cpp
	src/core/geometry/mesh_clip.cpp
	src/core/geometry/mesh_selection.cpp
	src/core/geometry/embedded_boundary.cpp
	src/core/geometry/local_surface_arrangement.cpp
	src/core/geometry/naca_theory.cpp
)
set_target_properties(libparacfd PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_compile_features(libparacfd PUBLIC cxx_std_20)
if(PARACFD_VALIDATION_FP64)
	target_compile_definitions(libparacfd PUBLIC PARACFD_REAL_DOUBLE=1)
endif()
# cudart_static is built against the static CRT; the rest of the build is /MD.
# Prefer the dynamic CRT (MSVCRT) everywhere to silence the benign LNK4098.
if(MSVC)
	target_link_options(libparacfd INTERFACE /NODEFAULTLIB:LIBCMT)
endif()
# Public: our own headers (include as "core/...") and the vendored single header.
target_include_directories(libparacfd PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty
	${CUDAToolkit_INCLUDE_DIRS}
)

# --- paracfd_geometry (STEP import via OpenCascade; isolates OCC from the core) -
# Qt-free static lib that mirrors how paracfd_gui_cuda splits out the GL/CUDA interop:
# it is the ONLY target that links OpenCascade, so libparacfd and every physics gate
# test stay OCC-free. Only paracfd-gui and the geometry unit test link it. STEP files
# arrive in millimetres and are scaled to metres inside step_import.cpp (PLAN §3).
set(OPENCASCADE_DIR "C:/OpenCASCADE-8.0/build2" CACHE PATH "OpenCascade 8.0 install root")
set(OCC_INCLUDE_DIR "${OPENCASCADE_DIR}/inc")
set(OCC_LIB_DIR "${OPENCASCADE_DIR}/win64/vc14/lib")
set(OCC_BIN_DIR "${OPENCASCADE_DIR}/win64/vc14/bin")

# Minimal toolkits needed to READ + TRIANGULATE a STEP shape (PLAN §3). These are the
# .lib import libs we link; being PRIVATE on a STATIC lib, CMake still propagates them
# (as $<LINK_ONLY>) to paracfd-gui / the geometry test so their final link resolves.
set(OCC_LINK_TOOLKITS
	TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKMesh TKPrim
	TKBO TKBool TKXSBase TKDESTEP)

# Full transitive runtime DLL closure of the link set (computed via dumpbin /dependents:
# TKDESTEP.dll pulls in TKDE/TKXCAF/TKLCAF/TKShHealing → TKV3d/TKService → freetype, etc.).
# These + the 3rdparty DLLs must sit next to any exe that loads paracfd_geometry so it stays
# double-clickable with NO OpenCascade on PATH (see CLAUDE.md).
set(OCC_RUNTIME_TOOLKITS
	TKBO TKBRep TKCAF TKCDF TKDE TKDESTEP TKernel TKG2d TKG3d TKGeomAlgo TKGeomBase
	TKHLR TKLCAF TKMath TKMesh TKPrim TKService TKShHealing TKTopAlgo TKV3d TKVCAF TKXCAF TKXSBase)
set(OCC_3RDPARTY_DLLS freetype brotlicommon brotlidec bz2 libpng16 zlib1)

# Copy the OpenCascade runtime DLLs next to a target's exe (POST_BUILD, copy_if_different so
# it is a no-op once deployed). Called for both paracfd-gui and the geometry test exe.
function(paracfd_deploy_occ_dlls target)
	if(NOT WIN32)
		return()
	endif()
	set(_paths)
	foreach(_tk ${OCC_RUNTIME_TOOLKITS} ${OCC_3RDPARTY_DLLS})
		set(_p "${OCC_BIN_DIR}/${_tk}.dll")
		if(EXISTS "${_p}")
			list(APPEND _paths "${_p}")
		else()
			message(WARNING "paracfd_deploy_occ_dlls(${target}): OCC DLL not found: ${_p}")
		endif()
	endforeach()
	add_custom_command(TARGET ${target} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_paths} "$<TARGET_FILE_DIR:${target}>"
		COMMENT "Deploying OpenCascade runtime DLLs next to ${target} (double-clickable, no OCC on PATH)"
		VERBATIM)
endfunction()

add_library(paracfd_geometry STATIC
	src/core/geometry/step_import.cpp
	src/core/geometry/step_face_uv_projection.cpp
	src/core/geometry/exact_cell_decomposition.cpp
	src/core/geometry/uv_constraint_mesh.cpp
	src/core/geometry/occt_surface_normalization.cpp
	src/core/geometry/naca_step.cpp)
set_target_properties(paracfd_geometry PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(paracfd_geometry PUBLIC cxx_std_20)
# OCC headers are PRIVATE — consumers include only the OCC-free step_import.h.
target_include_directories(paracfd_geometry
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
	PRIVATE ${OCC_INCLUDE_DIR})
# OCCT 8.0 moved some typedefs into a Deprecated/NCollectionAliases dir; include it +
# silence the deprecation notices (matches the OCC 8.0 usage in cobod-slicer).
set(OCC_DEPRECATED_DIR "${OPENCASCADE_DIR}/../OCCT/src/Deprecated/NCollectionAliases")
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(paracfd_geometry PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(paracfd_geometry PRIVATE OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
foreach(_tk ${OCC_LINK_TOOLKITS})
	if(NOT EXISTS "${OCC_LIB_DIR}/${_tk}.lib")
		message(FATAL_ERROR "paracfd_geometry: OpenCascade import lib not found: ${OCC_LIB_DIR}/${_tk}.lib")
	endif()
	target_link_libraries(paracfd_geometry PRIVATE "${OCC_LIB_DIR}/${_tk}.lib")
endforeach()

add_executable(paraglider_probe tools/paraglider_probe.cpp)
target_link_libraries(paraglider_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(paraglider_probe)

# Reproducible analytic CAD for aerodynamic validation, independent of imported
# paraglider topology and normals. The output is a regular STEP file consumed by
# exactly the same OpenCascade/EB pipeline as a wing from another application.
add_executable(naca_step_generator tools/naca_step_generator.cpp)
target_link_libraries(naca_step_generator PRIVATE paracfd_geometry)
paracfd_deploy_occ_dlls(naca_step_generator)

# Full imported-wing convergence/performance diagnostic. This is intentionally not a
# CTest because realistic STEP preprocessing and multi-step GPU runs are comparatively slow.
add_executable(paraglider_case_probe tools/paraglider_case_probe.cpp)
target_link_libraries(paraglider_case_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(paraglider_case_probe)

# Deliberately explicit/long-running validation executable rather than a default
# CTest: four profiles x three incidences plus half-cell phase checks take tens of
# minutes at the accepted 256-cells-per-chord finest spacing.
add_executable(naca_validation_probe tools/naca_validation_probe.cpp)
target_link_libraries(naca_validation_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(naca_validation_probe)

# --- parity_probe: immutable CPU/GPU oracle for reusable FP64 MAC kernels -----
enable_testing()

# Verify that OpenCascade import preserves the topology/precision contract needed by
# the zero-thickness EB classifier across tessellation settings.  This intentionally
# links only the isolated geometry target: no OCC types leak into the CFD core.
add_executable(step_import_provenance_probe
	tools/step_import_provenance_probe.cpp
	src/core/geometry/mesh_clip.cpp)
target_link_libraries(step_import_provenance_probe PRIVATE paracfd_geometry)
paracfd_deploy_occ_dlls(step_import_provenance_probe)
add_test(NAME step_import_provenance COMMAND step_import_provenance_probe)
set_tests_properties(step_import_provenance PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# OCC-free face-local constrained triangulation. The probe exercises the production
# geometry component without requiring a STEP fixture or OpenCascade type in its API.
add_executable(uv_constraint_mesh_probe tools/uv_constraint_mesh_probe.cpp)
target_link_libraries(uv_constraint_mesh_probe PRIVATE paracfd_geometry)
add_test(NAME uv_constraint_mesh COMMAND uv_constraint_mesh_probe)

# Direct OCCT regression for deterministic trim-valid contact UV charts. It
# manufactures shifted/full-period cylinders to exercise exact boundary
# pcurves, periodic lifts, 3D round trips, ambiguity and seam discontinuity.
add_executable(step_face_uv_projection_probe
	tools/step_face_uv_projection_probe.cpp)
target_compile_features(step_face_uv_projection_probe PRIVATE cxx_std_20)
target_include_directories(step_face_uv_projection_probe PRIVATE ${OCC_INCLUDE_DIR})
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(step_face_uv_projection_probe PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(step_face_uv_projection_probe PRIVATE
	OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
target_link_libraries(step_face_uv_projection_probe PRIVATE paracfd_geometry)
paracfd_deploy_occ_dlls(step_face_uv_projection_probe)
add_test(NAME step_face_uv_projection COMMAND step_face_uv_projection_probe)

# These topology checks intentionally use local CAD exports that are not repository
# fixtures. Keep the clean-checkout test graph independent of them while allowing the
# corrected and known-bad A/B models to be checked explicitly during development.
option(PARACFD_ENABLE_LOCAL_STEP_TESTS
	"Register topology tests that require local user STEP fixtures" OFF)
if(PARACFD_ENABLE_LOCAL_STEP_TESTS)
	foreach(_local_step_fixture
		"${CMAKE_SOURCE_DIR}/Test-Data/Two-Cells.stp"
		"${CMAKE_SOURCE_DIR}/Test-Data/Two-Cells-fixed.stp"
		"${CMAKE_SOURCE_DIR}/Test-Data/MiniripFixed.step")
		if(NOT EXISTS "${_local_step_fixture}")
			message(FATAL_ERROR
				"PARACFD_ENABLE_LOCAL_STEP_TESTS requires ${_local_step_fixture}")
		endif()
	endforeach()
	add_test(NAME step_topology_known_bad COMMAND step_import_provenance_probe
		Test-Data/Two-Cells.stp --quality-only --expect-components 3
		--expect-warnings 2 --expect-contacts 47 --expect-source-edges 72
		--expect-uv-failures 6)
	set_tests_properties(step_topology_known_bad PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	add_test(NAME step_topology_fixed_small COMMAND step_import_provenance_probe
		Test-Data/Two-Cells-fixed.stp --quality-only --expect-components 1
		--expect-warnings 0 --expect-contacts 51 --expect-source-edges 76
		--expect-shared-nodes-min 1 --expect-uv-failures 6)
	set_tests_properties(step_topology_fixed_small PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	add_test(NAME step_topology_fixed_full COMMAND step_import_provenance_probe
		Test-Data/MiniripFixed.step --quality-only --expect-components 1
		--expect-warnings 0 --expect-contacts 235 --expect-source-edges 412
		--expect-shared-nodes-min 1 --expect-partial-overlaps 2
		--expect-uv-failures 34)
	set_tests_properties(step_topology_fixed_full PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

add_executable(mesh_selection_probe
	tools/mesh_selection_probe.cpp
	src/core/geometry/mesh_selection.cpp)
target_compile_features(mesh_selection_probe PRIVATE cxx_std_20)
target_include_directories(mesh_selection_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME mesh_selection COMMAND mesh_selection_probe)

add_executable(topology_vertex_probe tools/topology_vertex_probe.cpp)
target_link_libraries(topology_vertex_probe PRIVATE libparacfd)
add_test(NAME topology_vertex_ids COMMAND topology_vertex_probe)

# Lightweight imported-CAD half-edge topology diagnostic. It deliberately stops at
# STEP/BVH preprocessing so classifier failures can be counted without building EB/AMR.
add_executable(cad_edge_classifier_probe tools/cad_edge_classifier_probe.cpp)
target_link_libraries(cad_edge_classifier_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(cad_edge_classifier_probe)

# Standalone diagnostic for evaluating OCCT BRep face imprinting before tessellation.
# It deliberately does not link the CFD/embedded-boundary core and is not a production path.
add_executable(occt_topology_normalization_probe tools/occt_topology_normalization_probe.cpp)
target_compile_features(occt_topology_normalization_probe PRIVATE cxx_std_20)
target_include_directories(occt_topology_normalization_probe PRIVATE ${OCC_INCLUDE_DIR})
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(occt_topology_normalization_probe PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(occt_topology_normalization_probe PRIVATE OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
target_link_libraries(occt_topology_normalization_probe PRIVATE paracfd_geometry
	"${OCC_LIB_DIR}/TKBO.lib" "${OCC_LIB_DIR}/TKShHealing.lib")
paracfd_deploy_occ_dlls(occt_topology_normalization_probe)

# Transactional exact-CAD contact experiment.  Unlike the broad normalization probe,
# this measures two-face and one-cell imprints around certified rib/skin contacts and
# never sews or inflates tolerances.  It is diagnostic only, not a CTest/production path.
add_executable(occt_contact_imprint_probe tools/occt_contact_imprint_probe.cpp)
target_compile_features(occt_contact_imprint_probe PRIVATE cxx_std_20)
target_include_directories(occt_contact_imprint_probe PRIVATE ${OCC_INCLUDE_DIR})
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(occt_contact_imprint_probe PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(occt_contact_imprint_probe PRIVATE OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
target_link_libraries(occt_contact_imprint_probe PRIVATE paracfd_geometry
	"${OCC_LIB_DIR}/TKBO.lib" "${OCC_LIB_DIR}/TKShHealing.lib")
paracfd_deploy_occ_dlls(occt_contact_imprint_probe)

# Guarded face-local constrained-meshing spike.  It injects a manufactured interior
# T-junction by default and can audit every curved interior-contact face in a STEP
# model.  Every requested segment must be a two-sided edge in both the transient mesh
# and final Poly_Triangulation; real-mode success requires every candidate to pass.
add_executable(occt_constrained_mesh_probe tools/occt_constrained_mesh_probe.cpp)
target_compile_features(occt_constrained_mesh_probe PRIVATE cxx_std_20)
target_include_directories(occt_constrained_mesh_probe PRIVATE ${OCC_INCLUDE_DIR})
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(occt_constrained_mesh_probe PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(occt_constrained_mesh_probe PRIVATE
	OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
target_link_libraries(occt_constrained_mesh_probe PRIVATE paracfd_geometry)
paracfd_deploy_occ_dlls(occt_constrained_mesh_probe)
add_test(NAME occt_constrained_mesh COMMAND occt_constrained_mesh_probe)

# Production planar-overlap normalizer probe. The helper retains one-to-many source
# provenance and proves coverage/load-area conservation before accepting each cluster;
# STEP-import integration is a separate, subsequent milestone.
add_executable(occt_overlap_normalization_probe
	tools/occt_overlap_normalization_probe.cpp)
target_compile_features(occt_overlap_normalization_probe PRIVATE cxx_std_20)
target_include_directories(occt_overlap_normalization_probe PRIVATE ${OCC_INCLUDE_DIR})
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(occt_overlap_normalization_probe PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(occt_overlap_normalization_probe PRIVATE
	OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
target_link_libraries(occt_overlap_normalization_probe PRIVATE paracfd_geometry
	"${OCC_LIB_DIR}/TKBO.lib")
paracfd_deploy_occ_dlls(occt_overlap_normalization_probe)
add_test(NAME occt_overlap_normalization_manufactured
	COMMAND occt_overlap_normalization_probe --manufactured-only)
if(PARACFD_ENABLE_LOCAL_STEP_TESTS)
	add_test(NAME occt_overlap_normalization_fixed_small
		COMMAND occt_overlap_normalization_probe --step Test-Data/Two-Cells-fixed.stp
			--expect-clusters 1 --expect-overlap-pairs 1)
	set_tests_properties(occt_overlap_normalization_fixed_small PROPERTIES
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	add_test(NAME occt_overlap_normalization_fixed_full
		COMMAND occt_overlap_normalization_probe --step Test-Data/MiniripFixed.step
			--expect-clusters 2 --expect-overlap-pairs 2)
	set_tests_properties(occt_overlap_normalization_fixed_full PROPERTIES
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

# OCCT-backed reference decomposition of finite, conforming zero-thickness sheets inside
# one Cartesian cell. This is the correctness gate for replacing the inference-heavy EB
# topology builder; it stays in the OCC-isolated preprocessing target and never runs per step.
add_executable(occt_cell_split_probe tools/occt_cell_split_probe.cpp)
target_link_libraries(occt_cell_split_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(occt_cell_split_probe)
add_test(NAME exact_cell_decomposition COMMAND occt_cell_split_probe)

# Fail-closed transaction gate for exact-cell failures at AMR ownership seams and
# for shared faces whose second participant rejects after an earlier face committed.
add_executable(amr_eb_transaction_probe tools/amr_eb_transaction_probe.cpp)
target_link_libraries(amr_eb_transaction_probe PRIVATE libparacfd)
add_test(NAME amr_eb_transactions COMMAND amr_eb_transaction_probe)

# A three-node aperture fan can contain a real sub-threshold fluid connector for
# which either local merge would reverse the surviving pressure distance. Exercise
# the production atlas stabilizer through its deterministic exact-cell test seam.
add_executable(amr_small_root_probe tools/amr_small_root_probe.cpp)
target_link_libraries(amr_small_root_probe PRIVATE libparacfd)
add_test(NAME amr_small_root COMMAND amr_small_root_probe)

# The strict composite pressure builder rejects an unsupported compact WLS
# correction. Only the explicitly labelled qualitative preview may fall back to
# the conservative orthogonal A/d term, and that fallback remains observable.
add_executable(amr_pressure_preview_probe tools/amr_pressure_preview_probe.cpp)
target_link_libraries(amr_pressure_preview_probe PRIVATE libparacfd)
add_test(NAME amr_pressure_preview_fallback COMMAND amr_pressure_preview_probe)

# Freeze the currently shipped qualitative-preview EB behavior before its legacy
# inference path is quarantined/replaced. This is intentionally a separate test from
# the strict geometry oracle and must never be interpreted as CAD certification.
add_executable(preview_topology_probe tools/preview_topology_probe.cpp)
target_link_libraries(preview_topology_probe PRIVATE libparacfd)
add_test(NAME preview_topology COMMAND preview_topology_probe)

add_executable(parity_probe tools/parity_probe.cu)
set_target_properties(parity_probe PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF)
target_link_libraries(parity_probe PRIVATE libparacfd)
add_test(NAME parity COMMAND parity_probe)
set_tests_properties(parity PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# Pure CPU oracle for the bounded single-run/AoA-sweep convergence policy. Keeping
# this independent of Qt, CUDA, and the geometry library makes exit ordering and
# final-window averaging deterministic and cheap to test.
add_executable(aero_convergence_probe tools/aero_convergence_probe.cpp)
target_include_directories(aero_convergence_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME aero_convergence COMMAND aero_convergence_probe)

# Deterministic CPU geometry/EB/AMR gates plus a tiny CUDA halo parity check. This is the
# zero-thickness membrane acceptance oracle; it does not use snapshot goldens.
add_executable(paraglider_geometry_probe tools/paraglider_geometry_probe.cpp)
target_link_libraries(paraglider_geometry_probe PRIVATE libparacfd paracfd_geometry)
paracfd_deploy_occ_dlls(paraglider_geometry_probe)
add_test(NAME paraglider_geometry COMMAND paraglider_geometry_probe)
set_tests_properties(paraglider_geometry PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# Imported analytic airfoil regression: the closed NACA extrusion must retain a
# separate internal fluid pressure component with no dead pressure mappings.
add_test(NAME naca_geometry COMMAND paraglider_probe --config configs/naca2412.json
	--expect-gauges 1 --expect-no-inactive-mappings)
set_tests_properties(naca_geometry PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

add_executable(paraglider_gpu_probe tools/paraglider_gpu_probe.cpp)
target_link_libraries(paraglider_gpu_probe PRIVATE libparacfd)
add_test(NAME paraglider_gpu COMMAND paraglider_gpu_probe)

add_executable(paraglider_flow_probe tools/paraglider_flow_probe.cpp)
target_link_libraries(paraglider_flow_probe PRIVATE libparacfd)
add_test(NAME paraglider_flow COMMAND paraglider_flow_probe)
set_tests_properties(paraglider_flow PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

add_executable(local_surface_arrangement_probe tools/local_surface_arrangement_probe.cpp)
target_link_libraries(local_surface_arrangement_probe PRIVATE libparacfd)
add_test(NAME local_surface_arrangement COMMAND local_surface_arrangement_probe)

# --- paracfd_gui_cuda (Qt-free slice sampler/colour reference) ----------------
add_library(paracfd_gui_cuda STATIC
	src/gui/slice_field.cu   # shared sampler/colourmap CPU reference (+ CUDA validation kernel)
)
set_target_properties(paracfd_gui_cuda PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(paracfd_gui_cuda PUBLIC libparacfd)

add_executable(visualization_probe tools/visualization_probe.cpp)
target_link_libraries(visualization_probe PRIVATE paracfd_gui_cuda)
add_test(NAME visualization_fields COMMAND visualization_probe)

# --- paracfd-gui (Qt6 + GL paraglider viewer) ---------------------------------
# Default ON. Turn OFF to build
# a Qt-less stub (main_stub.cpp) on hosts without Qt.
option(PARACFD_ENABLE_QT "Build the Qt6 ParaCFD paraglider application" ON)
add_executable(paracfd-gui)
if(PARACFD_ENABLE_QT)
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
	message(STATUS "paracfd-gui: Qt6 from ${_qt_found}")
	target_sources(paracfd-gui PRIVATE
		src/gui/paraglider_main.cpp
		src/gui/paraglider_window.cpp
		src/gui/slice_viewer.cpp
		src/gui/paraglider_sim_worker.cpp
		src/gui/flow_particles.cpp
		src/gui/flow_tracers.cpp
	)
	set_target_properties(paracfd-gui PROPERTIES AUTOMOC ON)
	target_link_libraries(paracfd-gui PRIVATE
		libparacfd paracfd_gui_cuda paracfd_geometry
		Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
	)
	target_compile_definitions(paracfd-gui PRIVATE PARACFD_HAVE_QT=1)

	# paracfd_geometry links OpenCascade DLLs -> deploy them next to paracfd-gui.exe so the
	# exe stays double-clickable (windeployqt only handles Qt, not OCC).
	paracfd_deploy_occ_dlls(paracfd-gui)

	# --- Deploy the Qt runtime next to the exe so paracfd-gui.exe is double-clickable
	# (no Qt on PATH needed) and stays so after a clean rebuild. windeployqt copies the
	# Qt6 DLLs + the platforms/qwindows.dll plugin into the output dir. Runs after each
	# link; it is a no-op when the deployed files are already up to date.
	find_program(WINDEPLOYQT_EXECUTABLE windeployqt HINTS "${_qt_found}/bin")
	if(WINDEPLOYQT_EXECUTABLE)
		add_custom_command(TARGET paracfd-gui POST_BUILD
			COMMAND "${WINDEPLOYQT_EXECUTABLE}"
				"$<IF:$<CONFIG:Debug>,--debug,--release>"
				--no-translations --no-system-d3d-compiler "$<TARGET_FILE:paracfd-gui>"
			COMMENT "windeployqt: deploying Qt runtime next to paracfd-gui.exe (double-clickable)"
			VERBATIM)
		# windeployqt only ships the windows platform plugin; also deploy qoffscreen so
		# the headless --offscreen smoke path (STEP-load CI) works with no Qt on PATH.
		set(_qoffscreen "${_qt_found}/plugins/platforms/$<IF:$<CONFIG:Debug>,qoffscreend,qoffscreen>.dll")
		add_custom_command(TARGET paracfd-gui POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:paracfd-gui>/platforms"
			COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_qoffscreen}" "$<TARGET_FILE_DIR:paracfd-gui>/platforms"
			COMMENT "Deploying Qt offscreen platform plugin next to paracfd-gui.exe"
			VERBATIM)
	else()
		message(WARNING "paracfd-gui: windeployqt not found under ${_qt_found}/bin; "
			"the exe will need the Qt bin dir on PATH to run (see CLAUDE.md).")
	endif()

else()
	target_sources(paracfd-gui PRIVATE src/gui/main_stub.cpp)
	target_link_libraries(paracfd-gui PRIVATE libparacfd)
endif()
