cmake_minimum_required(VERSION 3.24)

# WindCFD — 3D morphodynamic scour-protection simulator (C++/CUDA + Qt6).
# Milestone M0: scaffold. Three targets (libwindcfd, scour CLI, windcfd-gui stub),
# GoogleTest via FetchContent, ctest wiring, JSON config loader, VTI writer.
project(WindCFD LANGUAGES CXX CUDA)

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
# path explicitly; exported PUBLIC from libwindcfd so all consumers pick it up.
find_package(CUDAToolkit REQUIRED)

# --- libwindcfd (static core, no Qt) -----------------------------------------
add_library(libwindcfd STATIC
	src/core/config.cpp
	src/core/vti_writer.cpp
	src/core/vti_reader.cpp
	# Wind loads — Cp + force/moment integration over the voxelized building.
	src/core/windloads.cpp
	src/core/cuda_probe.cu
	# M1 — fluid core (MAC + MacCormack + Smagorinsky + MGPCG)
	src/core/fluid/advect.cu
	src/core/fluid/turbulence.cu
	src/core/fluid/project.cu
	src/core/fluid/mgpcg.cu
	src/core/fluid/stam_fluid_core.cu
	src/core/fluid/cavity.cpp
	# M2 — open channel + obstacles (inlet/Orlanski/lid, voxelized obstacle mask)
	src/core/fluid/channel_ops.cu
	src/core/fluid/channel_pressure.cu
	src/core/fluid/channel_core.cu
	src/core/fluid/channel_porous.cu
	src/core/fluid/channel_empty.cpp
	# M3 — bed-shear wall function, periodic channel + body force, SEM inlet, precursor
	src/core/fluid/bedshear.cu
	src/core/fluid/channel_periodic.cu
	src/core/fluid/plane_ops.cu
	src/core/fluid/sem_inlet.cu
	src/core/fluid/channel_sem.cpp
	src/core/fluid/precursor.cpp
	# Geometry — watertight ray-parity voxelizer (STEP TriMesh → solid obstacle mask). OCC-free
	# (needs only the mesh); the shared model placement + TriMesh headers live in core/geometry/.
	src/core/geometry/voxelize.cpp
	# Building — 3D-print centerline surface -> thickened wall + flat roof solid mask (distance-field).
	src/core/geometry/building.cpp
)
set_target_properties(libwindcfd PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_compile_features(libwindcfd PUBLIC cxx_std_20)
# cudart_static is built against the static CRT; the rest of the build is /MD.
# Prefer the dynamic CRT (MSVCRT) everywhere to silence the benign LNK4098.
if(MSVC)
	target_link_options(libwindcfd INTERFACE /NODEFAULTLIB:LIBCMT)
endif()
# Public: our own headers (include as "core/...") and the vendored single header.
target_include_directories(libwindcfd PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty
	${CUDAToolkit_INCLUDE_DIRS}
)

# --- windcfd_geometry (STEP import via OpenCascade; isolates OCC from the core) -
# Qt-free static lib that mirrors how windcfd_gui_cuda splits out the GL/CUDA interop:
# it is the ONLY target that links OpenCascade, so libwindcfd and every physics gate
# test stay OCC-free. Only windcfd-gui and the geometry unit test link it. STEP files
# arrive in millimetres and are scaled to metres inside step_import.cpp (PLAN §3).
set(OPENCASCADE_DIR "C:/OpenCASCADE-8.0/build2" CACHE PATH "OpenCascade 8.0 install root")
set(OCC_INCLUDE_DIR "${OPENCASCADE_DIR}/inc")
set(OCC_LIB_DIR "${OPENCASCADE_DIR}/win64/vc14/lib")
set(OCC_BIN_DIR "${OPENCASCADE_DIR}/win64/vc14/bin")

# Minimal toolkits needed to READ + TRIANGULATE a STEP shape (PLAN §3). These are the
# .lib import libs we link; being PRIVATE on a STATIC lib, CMake still propagates them
# (as $<LINK_ONLY>) to windcfd-gui / the geometry test so their final link resolves.
set(OCC_LINK_TOOLKITS
	TKernel TKMath TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep TKTopAlgo TKMesh TKXSBase TKDESTEP)

# Full transitive runtime DLL closure of the link set (computed via dumpbin /dependents:
# TKDESTEP.dll pulls in TKDE/TKXCAF/TKLCAF/TKShHealing → TKV3d/TKService → freetype, etc.).
# These + the 3rdparty DLLs must sit next to any exe that loads windcfd_geometry so it stays
# double-clickable with NO OpenCascade on PATH (see CLAUDE.md).
set(OCC_RUNTIME_TOOLKITS
	TKBO TKBRep TKCAF TKCDF TKDE TKDESTEP TKernel TKG2d TKG3d TKGeomAlgo TKGeomBase
	TKHLR TKLCAF TKMath TKMesh TKPrim TKService TKShHealing TKTopAlgo TKV3d TKVCAF TKXCAF TKXSBase)
set(OCC_3RDPARTY_DLLS freetype brotlicommon brotlidec bz2 libpng16 zlib1)

# Copy the OpenCascade runtime DLLs next to a target's exe (POST_BUILD, copy_if_different so
# it is a no-op once deployed). Called for both windcfd-gui and the geometry test exe.
function(windcfd_deploy_occ_dlls target)
	if(NOT WIN32)
		return()
	endif()
	set(_paths)
	foreach(_tk ${OCC_RUNTIME_TOOLKITS} ${OCC_3RDPARTY_DLLS})
		set(_p "${OCC_BIN_DIR}/${_tk}.dll")
		if(EXISTS "${_p}")
			list(APPEND _paths "${_p}")
		else()
			message(WARNING "windcfd_deploy_occ_dlls(${target}): OCC DLL not found: ${_p}")
		endif()
	endforeach()
	add_custom_command(TARGET ${target} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_paths} "$<TARGET_FILE_DIR:${target}>"
		COMMENT "Deploying OpenCascade runtime DLLs next to ${target} (double-clickable, no OCC on PATH)"
		VERBATIM)
endfunction()

add_library(windcfd_geometry STATIC src/core/geometry/step_import.cpp)
set_target_properties(windcfd_geometry PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(windcfd_geometry PUBLIC cxx_std_20)
# OCC headers are PRIVATE — consumers include only the OCC-free step_import.h.
target_include_directories(windcfd_geometry
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
	PRIVATE ${OCC_INCLUDE_DIR})
# OCCT 8.0 moved some typedefs into a Deprecated/NCollectionAliases dir; include it +
# silence the deprecation notices (matches the OCC 8.0 usage in cobod-slicer).
set(OCC_DEPRECATED_DIR "${OPENCASCADE_DIR}/../OCCT/src/Deprecated/NCollectionAliases")
if(EXISTS "${OCC_DEPRECATED_DIR}")
	target_include_directories(windcfd_geometry PRIVATE ${OCC_DEPRECATED_DIR})
endif()
target_compile_definitions(windcfd_geometry PRIVATE OCCT_NO_DEBUG OCCT_NO_DEPRECATED)
foreach(_tk ${OCC_LINK_TOOLKITS})
	if(NOT EXISTS "${OCC_LIB_DIR}/${_tk}.lib")
		message(FATAL_ERROR "windcfd_geometry: OpenCascade import lib not found: ${OCC_LIB_DIR}/${_tk}.lib")
	endif()
	target_link_libraries(windcfd_geometry PRIVATE "${OCC_LIB_DIR}/${_tk}.lib")
endforeach()

# --- building_probe (dev tool: centerline STEP -> building solid mask) -------
# Loads a centerline STEP via windcfd_geometry, sections + voxelizes a building. Links OCC
# (through windcfd_geometry) so the OCC DLLs are deployed next to it.
add_executable(building_probe tools/building_probe.cpp)
target_link_libraries(building_probe PRIVATE libwindcfd windcfd_geometry)
windcfd_deploy_occ_dlls(building_probe)

# --- parity_probe (Phase-0 migration oracle for the graded-structured-grid change) ---
# Restores the CLAUDE.md "1e-5 GPU-vs-CPU parity per kernel" invariant (no test target had
# existed — the CPU twins were defined but never run) AND adds the behaviour-preservation
# "grading collapses to uniform" golden gate. For every MAC-solver kernel in the Phase-A
# rewrite scope it checks (1) GPU kernel == CPU twin at rel. max-norm 1e-5 and (2) GPU output
# == a blessed snapshot. Bless once against the current uniform solver; every later run must
# still match, making "zero behaviour change" falsifiable. Pure CUDA + libwindcfd (no Qt/OCC/GL).
enable_testing()
add_executable(parity_probe tools/parity_probe.cu)
set_target_properties(parity_probe PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF)
target_link_libraries(parity_probe PRIVATE libwindcfd)
add_test(NAME parity COMMAND parity_probe)
set_tests_properties(parity PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# --- windcfd_gui_cuda (GUI CUDA half: slice colour kernel + CUDA-GL interop) ---
# Kept Qt-free so nvcc never sees Qt's host flags, and so the fast unit suite can
# link the slice kernel + its CPU reference without pulling in Qt/GL.
add_library(windcfd_gui_cuda STATIC
	src/gui/slice_field.cu   # pure CUDA sampler+colourmap (+ CPU reference)
	src/gui/slice_gl.cu      # cudaGraphicsGLRegisterBuffer interop (the only GL-header TU)
)
set_target_properties(windcfd_gui_cuda PROPERTIES
	CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
	CUDA_SEPARABLE_COMPILATION OFF
	POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(windcfd_gui_cuda PUBLIC libwindcfd)

# --- windcfd-gui (Qt6 + GL 4.3 slice viewer; G1) ------------------------------
# Default ON from G1 onward (Qt 6.11.1 present on this machine). Turn OFF to build
# a Qt-less stub (main_stub.cpp) on hosts without Qt.
option(WINDCFD_ENABLE_QT "Build windcfd-gui against Qt6 (G1 slice viewer)" ON)
add_executable(windcfd-gui)
if(WINDCFD_ENABLE_QT)
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
	message(STATUS "windcfd-gui: Qt6 from ${_qt_found}")
	target_sources(windcfd-gui PRIVATE
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
	set_target_properties(windcfd-gui PROPERTIES AUTOMOC ON)
	target_link_libraries(windcfd-gui PRIVATE
		libwindcfd windcfd_gui_cuda windcfd_geometry
		Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
	)
	target_compile_definitions(windcfd-gui PRIVATE WINDCFD_HAVE_QT=1)

	# windcfd_geometry links OpenCascade DLLs -> deploy them next to windcfd-gui.exe so the
	# exe stays double-clickable (windeployqt only handles Qt, not OCC).
	windcfd_deploy_occ_dlls(windcfd-gui)

	# --- Deploy the Qt runtime next to the exe so windcfd-gui.exe is double-clickable
	# (no Qt on PATH needed) and stays so after a clean rebuild. windeployqt copies the
	# Qt6 DLLs + the platforms/qwindows.dll plugin into the output dir. Runs after each
	# link; it is a no-op when the deployed files are already up to date.
	find_program(WINDEPLOYQT_EXECUTABLE windeployqt HINTS "${_qt_found}/bin")
	if(WINDEPLOYQT_EXECUTABLE)
		add_custom_command(TARGET windcfd-gui POST_BUILD
			COMMAND "${WINDEPLOYQT_EXECUTABLE}"
				"$<IF:$<CONFIG:Debug>,--debug,--release>"
				--no-translations --no-system-d3d-compiler "$<TARGET_FILE:windcfd-gui>"
			COMMENT "windeployqt: deploying Qt runtime next to windcfd-gui.exe (double-clickable)"
			VERBATIM)
		# windeployqt only ships the windows platform plugin; also deploy qoffscreen so
		# the headless --offscreen smoke path (STEP-load CI) works with no Qt on PATH.
		set(_qoffscreen "${_qt_found}/plugins/platforms/$<IF:$<CONFIG:Debug>,qoffscreend,qoffscreen>.dll")
		add_custom_command(TARGET windcfd-gui POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:windcfd-gui>/platforms"
			COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_qoffscreen}" "$<TARGET_FILE_DIR:windcfd-gui>/platforms"
			COMMENT "Deploying Qt offscreen platform plugin next to windcfd-gui.exe"
			VERBATIM)
	else()
		message(WARNING "windcfd-gui: windeployqt not found under ${_qt_found}/bin; "
			"the exe will need the Qt bin dir on PATH to run (see CLAUDE.md).")
	endif()

else()
	target_sources(windcfd-gui PRIVATE src/gui/main_stub.cpp)
	target_link_libraries(windcfd-gui PRIVATE libwindcfd)
endif()
