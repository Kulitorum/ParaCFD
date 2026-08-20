// occt_conforming_import_assembly.h -- final exact-CAD import hand-off.
//
// This is the last OpenCascade-facing stage of the STEP geometry pipeline.  It
// converts the face-local, exactly conformed mesh into the OCCT-free TriMesh and
// StepCadContactGraph contracts used by BVH/EB preprocessing.  No geometric
// search, welding, healing, or re-triangulation is performed here.
#pragma once

#include "core/geometry/occt_conforming_mesh_builder.h"
#include "core/geometry/step_import.h"

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace paracfd::core
{
	// StepContactFaceUse currently represents one logical face use.  A periodic
	// seam has two distinct UV-chart occurrences but one logical face use with two
	// sectors.  Preserve those branches losslessly here; the public graph keeps the
	// deterministic first branch in StepContactFaceUse::sample_uv for compatibility.
	struct OcctConformingContactUvBranch
	{
		std::uint64_t contact_atom_id = TriMesh::kNoCadContactId;
		std::uint32_t source_face_id = std::numeric_limits<std::uint32_t>::max();
		std::uint32_t boundary_occurrence_id =
			OcctContactBoundaryOccurrence::no_occurrence_id;
		std::int8_t boundary_orientation = 0;
		StepContactUseKind kind = StepContactUseKind::trim_boundary;
		std::vector<std::array<double, 2>> sample_uv;
	};

	struct OcctConformingImportAssembly
	{
		TriMesh mesh;
		StepCadContactGraph contacts;
		// Includes fan-one/open-boundary atoms as well as public fan>=2 contacts.
		// This is an audit sidecar, not a declaration that fan-one is an aerodynamic
		// contact.  Periodic seams have one record per exact oriented occurrence.
		std::vector<OcctConformingContactUvBranch> contact_uv_branches;
	};

	// Convert OCCT model units (the ParaCFD STEP pipeline uses millimetres) to the
	// SI-metre hand-off.  The operation is transactional: `output` remains logically
	// unchanged on every validation or allocation failure and only `error` changes.
	//
	// Strong validation includes:
	//  * bit-identical positions for every shared topology ID;
	//  * globally unique IDs for every ordinary face-local vertex;
	//  * complete/consistent triangle-edge incidence and CAD boundary lineage;
	//  * contact half-edge incidence equal to every atom's exact fan degree; and
	//  * preservation of both UV branches of periodic seam occurrences.
	bool assemble_occt_conforming_import(const OcctContactTopology& topology,
		const OcctConformingMesh& conformed,
		OcctConformingImportAssembly& output, std::string& error);
}
