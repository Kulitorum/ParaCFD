// step_import.cpp — see step_import.h. The ONLY translation unit that touches OpenCascade.
//
// Pipeline (PLAN §3 recipe, reference impl cobod-slicer src/app/widgets/render/mesh.cpp:494):
//   STEPControl_Reader → ReadFile → TransferRoots → OneShape → BRepMesh_IncrementalMesh →
//   per-face BRep_Tool::Triangulation (apply TopLoc_Location) → 1-based→0-based indices.
//
// Three correctness-critical traps, all handled here:
//   (a) mm→m — OCC coordinates are millimetres; every node is multiplied by 0.001.
//   (b) winding — OCC stores triangles in the surface's natural (u,v) sense; the solid's
//       OUTWARD normal is that sense only when the face is FORWARD, so we swap two indices
//       when face.Orientation() == TopAbs_REVERSED. (This is the OPPOSITE branch from the
//       slicer, whose meshes are globally inverted — see PLAN §3 / CLAUDE.md.)
//   (c) smooth normals — per-vertex normals are the area-weighted average of incident face
//       normals (summing UN-normalised cross products naturally area-weights, then normalise).
//
// load_step_mesh() accumulates the whole shape; load_step_solids() reuses the SAME per-shape
// accumulation once PER SOLID for the PLAN G3 drop/settle convex-piece decomposition (each solid
// → one convex collision piece). Both share read_step_shape() (read+transfer+triangulate) and
// mesh_from_faces() (faces → TriMesh) so the two entry points stay byte-for-byte consistent.
#include "core/geometry/step_import.h"

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace windcfd::core
{
	namespace
	{
		void set_error(std::string* error, const std::string& msg)
		{
			if (error) *error = msg;
		}

		// Read + transfer a STEP file to its single shape and triangulate it in place (BRepMesh).
		// Returns false with *error set on any read/transfer/null-shape failure. deflection is in
		// OCC's native millimetres. May throw Standard_Failure (OCC) — callers wrap in try/catch.
		bool read_step_shape(const std::string& path, double deflection_mm, TopoDS_Shape& out, std::string* error)
		{
			STEPControl_Reader reader;
			IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());
			if (status != IFSelect_RetDone)
			{
				set_error(error, "STEPControl_Reader::ReadFile failed for '" + path + "'");
				return false;
			}

			Standard_Integer nroots = reader.TransferRoots();
			if (nroots <= 0)
			{
				set_error(error, "STEPControl_Reader::TransferRoots transferred no roots");
				return false;
			}

			out = reader.OneShape();
			if (out.IsNull())
			{
				set_error(error, "reader.OneShape() returned a null shape");
				return false;
			}

			// Triangulate (deflection is in OCC's native millimetres).
			if (deflection_mm <= 0.0) deflection_mm = 0.1;
			BRepMesh_IncrementalMesh mesher(out, deflection_mm);
			mesher.Perform();
			return true;
		}

		// Accumulate every triangulated FACE of `shape` into one TriMesh (metres, outward winding,
		// area-weighted per-vertex normals, bbox). `shape` must ALREADY be meshed (BRepMesh run on
		// it or an ancestor). Returns empty() if the shape contributes no triangles.
		TriMesh mesh_from_faces(const TopoDS_Shape& shape)
		{
			TriMesh mesh;

			// --- Accumulate faces into one flat mesh (metres) ----------------------
			constexpr double kMmToM = 0.001;
			for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
			{
				TopoDS_Face face = TopoDS::Face(exp.Current());
				TopLoc_Location loc;
				Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
				if (tri.IsNull() || tri->NbNodes() <= 0 || tri->NbTriangles() <= 0) continue;

				const bool reversed = (face.Orientation() == TopAbs_REVERSED);
				const std::uint32_t base = static_cast<std::uint32_t>(mesh.positions.size() / 3);
				const gp_Trsf trsf = loc.Transformation();

				for (Standard_Integer i = 1; i <= tri->NbNodes(); ++i)
				{
					gp_Pnt p = tri->Node(i);
					p.Transform(trsf);
					mesh.positions.push_back(static_cast<float>(p.X() * kMmToM));
					mesh.positions.push_back(static_cast<float>(p.Y() * kMmToM));
					mesh.positions.push_back(static_cast<float>(p.Z() * kMmToM));
				}

				for (Standard_Integer i = 1; i <= tri->NbTriangles(); ++i)
				{
					Standard_Integer n1, n2, n3;
					tri->Triangle(i).Get(n1, n2, n3);
					--n1; --n2; --n3; // OCC nodes are 1-based
					if (reversed) std::swap(n1, n2); // outward winding for a solid's face
					mesh.indices.push_back(base + static_cast<std::uint32_t>(n1));
					mesh.indices.push_back(base + static_cast<std::uint32_t>(n2));
					mesh.indices.push_back(base + static_cast<std::uint32_t>(n3));
				}
			}

			if (mesh.indices.empty()) return mesh; // empty() — caller reports the error/skips

			// --- Per-vertex normals: area-weighted average of incident face normals ----
			const std::size_t nverts = mesh.vertex_count();
			std::vector<double> acc(nverts * 3, 0.0);
			for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3)
			{
				const std::uint32_t i0 = mesh.indices[t + 0];
				const std::uint32_t i1 = mesh.indices[t + 1];
				const std::uint32_t i2 = mesh.indices[t + 2];
				const float* p0 = &mesh.positions[3 * i0];
				const float* p1 = &mesh.positions[3 * i1];
				const float* p2 = &mesh.positions[3 * i2];
				const double e1x = p1[0] - p0[0], e1y = p1[1] - p0[1], e1z = p1[2] - p0[2];
				const double e2x = p2[0] - p0[0], e2y = p2[1] - p0[1], e2z = p2[2] - p0[2];
				// Cross(e1,e2): magnitude == 2*area, so accumulating raw crosses area-weights.
				const double nx = e1y * e2z - e1z * e2y;
				const double ny = e1z * e2x - e1x * e2z;
				const double nz = e1x * e2y - e1y * e2x;
				for (std::uint32_t idx : { i0, i1, i2 })
				{
					acc[3 * idx + 0] += nx;
					acc[3 * idx + 1] += ny;
					acc[3 * idx + 2] += nz;
				}
			}

			mesh.normals.resize(nverts * 3);
			for (std::size_t v = 0; v < nverts; ++v)
			{
				double nx = acc[3 * v + 0], ny = acc[3 * v + 1], nz = acc[3 * v + 2];
				const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
				if (len > 1e-20) { nx /= len; ny /= len; nz /= len; }
				else { nx = 0.0; ny = 0.0; nz = 1.0; } // degenerate fallback: z-up
				mesh.normals[3 * v + 0] = static_cast<float>(nx);
				mesh.normals[3 * v + 1] = static_cast<float>(ny);
				mesh.normals[3 * v + 2] = static_cast<float>(nz);
			}

			// --- Bounding box (metres) --------------------------------------------------
			float lo[3] = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
			float hi[3] = { -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() };
			for (std::size_t v = 0; v < nverts; ++v)
				for (int c = 0; c < 3; ++c)
				{
					const float x = mesh.positions[3 * v + c];
					lo[c] = std::min(lo[c], x);
					hi[c] = std::max(hi[c], x);
				}
			mesh.bbox_min = { { lo[0], lo[1], lo[2] } };
			mesh.bbox_max = { { hi[0], hi[1], hi[2] } };

			return mesh;
		}
	}

	TriMesh load_step_mesh(const std::string& path, double deflection_mm, std::string* error)
	{
		try
		{
			TopoDS_Shape shape;
			if (!read_step_shape(path, deflection_mm, shape, error)) return TriMesh{};

			TriMesh mesh = mesh_from_faces(shape);
			if (mesh.empty())
			{
				set_error(error, "shape produced no triangulable faces");
				return TriMesh{};
			}
			if (error) error->clear();
			return mesh;
		}
		catch (const Standard_Failure& f)
		{
			set_error(error, std::string("OpenCascade exception: ") + f.GetMessageString());
			return TriMesh{};
		}
	}

	std::vector<TriMesh> load_step_solids(const std::string& path, double deflection_mm, std::string* error)
	{
		std::vector<TriMesh> solids;
		try
		{
			TopoDS_Shape shape;
			if (!read_step_shape(path, deflection_mm, shape, error)) return {};

			for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next())
			{
				TriMesh m = mesh_from_faces(exp.Current());
				if (!m.empty()) solids.push_back(std::move(m));
			}

			// Shell/face-only STEP (no TopAbs_SOLID): fall back to the whole shape as one piece.
			if (solids.empty())
			{
				TriMesh whole = mesh_from_faces(shape);
				if (!whole.empty()) solids.push_back(std::move(whole));
			}

			if (solids.empty())
			{
				set_error(error, "STEP produced no triangulable solids");
				return {};
			}
		}
		catch (const Standard_Failure& f)
		{
			set_error(error, std::string("OpenCascade exception: ") + f.GetMessageString());
			return {};
		}

		if (error) error->clear();
		return solids;
	}
}
