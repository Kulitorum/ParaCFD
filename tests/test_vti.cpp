// test_vti.cpp — VTI writer/reader round-trip + validation (the M0 gate core).
#include "core/vti_synthetic.h"
#include "core/vti_writer.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace scour::io;

namespace
{
	std::string temp_vti(const char* stem)
	{
		auto dir = std::filesystem::temp_directory_path();
		auto p = dir / (std::string("scour_m0_") + stem + ".vti");
		return p.string();
	}
}

// Write the canonical synthetic field, re-parse it, and require an exact match of
// geometry + every payload sample. This is the producer/consumer gate for M0.
TEST(VtiRoundTrip, WriteReadExact)
{
	const std::string path = temp_vti("roundtrip");
	VtiImage src = make_synthetic_image();

	std::string err;
	ASSERT_TRUE(write_vti(path, src, &err)) << err;

	VtiImage got;
	ASSERT_TRUE(read_vti(path, got, &err)) << err;

	EXPECT_EQ(got.nx, src.nx);
	EXPECT_EQ(got.ny, src.ny);
	EXPECT_EQ(got.nz, src.nz);
	EXPECT_DOUBLE_EQ(got.ox, src.ox);
	EXPECT_DOUBLE_EQ(got.oy, src.oy);
	EXPECT_DOUBLE_EQ(got.oz, src.oz);
	EXPECT_DOUBLE_EQ(got.sx, src.sx);
	EXPECT_DOUBLE_EQ(got.sy, src.sy);
	EXPECT_DOUBLE_EQ(got.sz, src.sz);

	ASSERT_EQ(got.point_fields.size(), 1u);
	ASSERT_EQ(got.point_fields[0].name, "synthetic");
	ASSERT_EQ(got.point_fields[0].data.size(), src.num_points());

	// Raw appended float32 is byte-preserved: expect bit-exact equality, and also
	// verify against the analytic formula recomputed independently.
	const auto& a = src.point_fields[0].data;
	const auto& b = got.point_fields[0].data;
	std::size_t mismatches = 0;
	for (int k = 0; k < src.nz; ++k)
		for (int j = 0; j < src.ny; ++j)
			for (int i = 0; i < src.nx; ++i)
			{
				const std::size_t idx = static_cast<std::size_t>(i) + static_cast<std::size_t>(src.nx) * (static_cast<std::size_t>(j) + static_cast<std::size_t>(src.ny) * static_cast<std::size_t>(k));
				const float expect = static_cast<float>(synthetic_field_value(i, j, k, src));
				if (a[idx] != expect || b[idx] != a[idx])
					++mismatches;
			}
	EXPECT_EQ(mismatches, 0u);

	std::error_code ec;
	std::filesystem::remove(path, ec);
}

// The reader must reject a structurally invalid file rather than crash / succeed.
TEST(VtiRoundTrip, RejectsGarbage)
{
	const std::string path = temp_vti("garbage");
	{
		std::FILE* f = std::fopen(path.c_str(), "wb");
		ASSERT_NE(f, nullptr);
		const char junk[] = "not a vtk file at all";
		std::fwrite(junk, 1, sizeof(junk) - 1, f);
		std::fclose(f);
	}
	VtiImage got;
	std::string err;
	EXPECT_FALSE(read_vti(path, got, &err));
	EXPECT_FALSE(err.empty());

	std::error_code ec;
	std::filesystem::remove(path, ec);
}

// Field-size mismatch must be caught by the writer, not silently truncated.
TEST(VtiRoundTrip, WriterRejectsBadFieldSize)
{
	VtiImage img;
	img.nx = 4; img.ny = 4; img.nz = 4;
	img.sx = img.sy = img.sz = 0.1;
	VtiField f;
	f.name = "bad";
	f.data.resize(10); // != 64
	img.point_fields.push_back(f);
	std::string err;
	EXPECT_FALSE(write_vti(temp_vti("badsize"), img, &err));
	EXPECT_FALSE(err.empty());
}
