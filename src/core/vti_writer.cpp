// vti_writer.cpp — hand-rolled VTK ImageData writer (appended raw binary).
#include "core/vti_writer.h"

#include <cstring>
#include <fstream>
#include <sstream>

namespace windcfd::io
{
	namespace
	{
		// Append the little-endian raw bytes of a trivially-copyable value.
		template <typename T>
		void append_le(std::vector<char>& buf, const T& v)
		{
			// x86 is little-endian; copy the object representation directly.
			const char* p = reinterpret_cast<const char*>(&v);
			buf.insert(buf.end(), p, p + sizeof(T));
		}
	}

	bool write_vti(const std::string& path, const VtiImage& img, std::string* err)
	{
		auto fail = [&](const std::string& m) { if (err) *err = m; return false; };

		if (img.nx < 1 || img.ny < 1 || img.nz < 1)
			return fail("write_vti: dimensions must be >= 1");
		const std::size_t np = img.num_points();
		if (img.point_fields.empty())
			return fail("write_vti: no point fields to write");
		for (const auto& f : img.point_fields)
		{
			if (f.data.size() != np)
				return fail("write_vti: field '" + f.name + "' size " + std::to_string(f.data.size()) + " != num_points " + std::to_string(np));
		}

		// Compute appended-data offsets. Each array is [UInt64 nbytes][nbytes data].
		const std::uint64_t hdr = sizeof(std::uint64_t);            // 8-byte length header
		const std::uint64_t databytes = static_cast<std::uint64_t>(np) * sizeof(float);
		std::vector<std::uint64_t> offsets(img.point_fields.size());
		std::uint64_t cursor = 0;
		for (std::size_t i = 0; i < img.point_fields.size(); ++i)
		{
			offsets[i] = cursor;
			cursor += hdr + databytes;
		}

		// XML header.
		std::ostringstream xml;
		xml << "<?xml version=\"1.0\"?>\n";
		xml << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
		xml << "  <ImageData WholeExtent=\"0 " << (img.nx - 1) << " 0 " << (img.ny - 1) << " 0 " << (img.nz - 1) << "\""
			<< " Origin=\"" << img.ox << " " << img.oy << " " << img.oz << "\""
			<< " Spacing=\"" << img.sx << " " << img.sy << " " << img.sz << "\">\n";
		xml << "    <Piece Extent=\"0 " << (img.nx - 1) << " 0 " << (img.ny - 1) << " 0 " << (img.nz - 1) << "\">\n";
		xml << "      <PointData Scalars=\"" << img.point_fields.front().name << "\">\n";
		for (std::size_t i = 0; i < img.point_fields.size(); ++i)
		{
			xml << "        <DataArray type=\"Float32\" Name=\"" << img.point_fields[i].name
				<< "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offsets[i] << "\"/>\n";
		}
		xml << "      </PointData>\n";
		xml << "      <CellData/>\n";
		xml << "    </Piece>\n";
		xml << "  </ImageData>\n";
		xml << "  <AppendedData encoding=\"raw\">\n_";

		// Raw appended payload.
		std::vector<char> payload;
		payload.reserve(static_cast<std::size_t>(cursor));
		for (const auto& f : img.point_fields)
		{
			append_le<std::uint64_t>(payload, databytes);
			const char* dp = reinterpret_cast<const char*>(f.data.data());
			payload.insert(payload.end(), dp, dp + databytes);
		}

		std::ofstream out(path, std::ios::binary);
		if (!out)
			return fail("write_vti: cannot open output file: " + path);
		const std::string head = xml.str();
		out.write(head.data(), static_cast<std::streamsize>(head.size()));
		out.write(payload.data(), static_cast<std::streamsize>(payload.size()));
		const std::string tail = "\n  </AppendedData>\n</VTKFile>\n";
		out.write(tail.data(), static_cast<std::streamsize>(tail.size()));
		out.flush();
		if (!out)
			return fail("write_vti: write error on: " + path);
		return true;
	}
}
