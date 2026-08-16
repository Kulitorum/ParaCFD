// vti_reader.cpp — re-parse + validate a .vti written by write_vti.
//
// A deliberately narrow parser for our own appended-raw / UInt64 / Float32 /
// LittleEndian layout. It validates the header structure and, critically, checks
// each appended array's UInt64 byte-count header against the declared point count.
#include "core/vti_writer.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace paracfd::io
{
	namespace
	{
		// Read an attribute value: name="value". Searches [from,to) of s. Returns
		// false if not found.
		bool attr(const std::string& s, std::size_t from, std::size_t to, const char* name, std::string& out)
		{
			const std::string key = std::string(name) + "=\"";
			std::size_t p = s.find(key, from);
			if (p == std::string::npos || p >= to)
				return false;
			p += key.size();
			std::size_t q = s.find('"', p);
			if (q == std::string::npos || q > to)
				return false;
			out = s.substr(p, q - p);
			return true;
		}

		bool parse_ints(const std::string& v, std::vector<long long>& out)
		{
			std::istringstream is(v);
			long long x;
			out.clear();
			while (is >> x)
				out.push_back(x);
			return !out.empty();
		}

		bool parse_doubles(const std::string& v, std::vector<double>& out)
		{
			std::istringstream is(v);
			double x;
			out.clear();
			while (is >> x)
				out.push_back(x);
			return !out.empty();
		}
	}

	bool read_vti(const std::string& path, VtiImage& out, std::string* err)
	{
		auto fail = [&](const std::string& m) { if (err) *err = m; return false; };

		std::ifstream in(path, std::ios::binary);
		if (!in)
			return fail("read_vti: cannot open " + path);
		std::string buf((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		if (buf.empty())
			return fail("read_vti: empty file " + path);

		if (buf.find("type=\"ImageData\"") == std::string::npos)
			return fail("read_vti: not a VTK ImageData file");

		std::string sval;
		if (!attr(buf, 0, buf.size(), "header_type", sval) || sval != "UInt64")
			return fail("read_vti: expected header_type=UInt64");
		if (!attr(buf, 0, buf.size(), "byte_order", sval) || sval != "LittleEndian")
			return fail("read_vti: expected byte_order=LittleEndian");

		// WholeExtent -> point dims.
		if (!attr(buf, 0, buf.size(), "WholeExtent", sval))
			return fail("read_vti: missing WholeExtent");
		std::vector<long long> ext;
		if (!parse_ints(sval, ext) || ext.size() != 6)
			return fail("read_vti: bad WholeExtent");
		out.nx = static_cast<int>(ext[1] - ext[0] + 1);
		out.ny = static_cast<int>(ext[3] - ext[2] + 1);
		out.nz = static_cast<int>(ext[5] - ext[4] + 1);
		if (out.nx < 1 || out.ny < 1 || out.nz < 1)
			return fail("read_vti: non-positive dimensions");

		std::vector<double> od, sd;
		if (!attr(buf, 0, buf.size(), "Origin", sval) || !parse_doubles(sval, od) || od.size() != 3)
			return fail("read_vti: bad Origin");
		if (!attr(buf, 0, buf.size(), "Spacing", sval) || !parse_doubles(sval, sd) || sd.size() != 3)
			return fail("read_vti: bad Spacing");
		out.ox = od[0]; out.oy = od[1]; out.oz = od[2];
		out.sx = sd[0]; out.sy = sd[1]; out.sz = sd[2];

		const std::size_t np = out.num_points();

		// Locate the appended raw blob: the byte just past the '_' after <AppendedData>.
		std::size_t ap = buf.find("<AppendedData");
		if (ap == std::string::npos)
			return fail("read_vti: missing AppendedData");
		std::size_t us = buf.find('_', ap);
		if (us == std::string::npos)
			return fail("read_vti: missing appended-data underscore");
		const std::size_t base = us + 1;

		// Collect DataArray tags (Name + offset), stop at </PointData>.
		std::size_t pd_end = buf.find("</PointData>");
		if (pd_end == std::string::npos)
			return fail("read_vti: missing PointData");

		out.point_fields.clear();
		std::size_t scan = buf.find("<PointData");
		if (scan == std::string::npos)
			return fail("read_vti: missing PointData open tag");
		while (true)
		{
			std::size_t da = buf.find("<DataArray", scan);
			if (da == std::string::npos || da > pd_end)
				break;
			std::size_t tag_end = buf.find('>', da);
			if (tag_end == std::string::npos)
				return fail("read_vti: unterminated DataArray tag");

			std::string name, type, offs;
			if (!attr(buf, da, tag_end, "Name", name))
				return fail("read_vti: DataArray missing Name");
			if (!attr(buf, da, tag_end, "type", type) || type != "Float32")
				return fail("read_vti: DataArray '" + name + "' must be Float32");
			if (!attr(buf, da, tag_end, "offset", offs))
				return fail("read_vti: DataArray '" + name + "' missing offset");
			const std::size_t offset = static_cast<std::size_t>(std::stoull(offs));

			// Read the UInt64 byte-count header, then the raw floats.
			const std::size_t hpos = base + offset;
			if (hpos + sizeof(std::uint64_t) > buf.size())
				return fail("read_vti: array '" + name + "' header past EOF");
			std::uint64_t nbytes = 0;
			std::memcpy(&nbytes, buf.data() + hpos, sizeof(nbytes));
			if (nbytes != static_cast<std::uint64_t>(np) * sizeof(float))
				return fail("read_vti: array '" + name + "' byte-count " + std::to_string(nbytes) + " != expected " + std::to_string(np * sizeof(float)));
			const std::size_t dpos = hpos + sizeof(std::uint64_t);
			if (dpos + nbytes > buf.size())
				return fail("read_vti: array '" + name + "' data past EOF");

			VtiField f;
			f.name = name;
			f.data.resize(np);
			std::memcpy(f.data.data(), buf.data() + dpos, static_cast<std::size_t>(nbytes));
			out.point_fields.push_back(std::move(f));

			scan = tag_end + 1;
		}

		if (out.point_fields.empty())
			return fail("read_vti: no DataArray found");
		return true;
	}
}
