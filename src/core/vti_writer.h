// vti_writer.h — hand-rolled VTK ImageData (.vti) writer + reader.
//
// Format: XML header with DataArrays declared format="appended", followed by a
// single raw appended-data blob. header_type="UInt64", byte_order="LittleEndian"
// (x86). Each appended array is [UInt64 byteCount][raw little-endian Float32 ...].
// Scalar point fields only (all M0 needs); this is deliberately not a full VTK
// implementation — it is a self-contained round-trippable writer + validator.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace paracfd::io
{
	struct VtiField
	{
		std::string name;
		std::vector<float> data; // length must equal VtiImage::num_points()
	};

	// A uniform (ImageData) grid with scalar point data.
	// nx/ny/nz are the number of POINTS along each axis (WholeExtent 0..n-1).
	struct VtiImage
	{
		int nx = 0, ny = 0, nz = 0;      // points per axis
		double ox = 0, oy = 0, oz = 0;   // origin (m)
		double sx = 0, sy = 0, sz = 0;   // spacing (m)
		std::vector<VtiField> point_fields;

		std::size_t num_points() const
		{
			return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
		}
	};

	// Write img to path. Returns false and fills *err on failure. Point-data
	// ordering is x-fastest then y then z: idx = i + nx*(j + ny*k).
	bool write_vti(const std::string& path, const VtiImage& img, std::string* err = nullptr);

	// Parse a .vti written by write_vti (also validates header consistency and the
	// per-array UInt64 byte-count against the declared point count). Returns false
	// and fills *err on any structural mismatch.
	bool read_vti(const std::string& path, VtiImage& out, std::string* err = nullptr);
}
