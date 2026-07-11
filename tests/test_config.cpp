// test_config.cpp — JSON config loader unit tests.
#include "core/config.h"

#include <gtest/gtest.h>

#include <stdexcept>

using scour::core::Config;

TEST(Config, DefaultsMatchResearchSpec)
{
	// Empty object -> canonical RESEARCH.md §1 defaults.
	Config c = Config::from_json_string("{}");
	EXPECT_DOUBLE_EQ(c.rho, 1027.0);
	EXPECT_DOUBLE_EQ(c.nu, 1.36e-6);
	EXPECT_DOUBLE_EQ(c.rho_s, 2650.0);
	EXPECT_DOUBLE_EQ(c.porosity, 0.36);
	EXPECT_DOUBLE_EQ(c.phi_repose_deg, 32.0);
	EXPECT_DOUBLE_EQ(c.domain_x, 10.0);
	EXPECT_DOUBLE_EQ(c.domain_z, 5.0);
	EXPECT_DOUBLE_EQ(c.voxel_h, 0.05);
}

TEST(Config, DerivedGridCounts)
{
	Config c = Config::from_json_string(R"({"domain_x":10.0,"domain_y":10.0,"domain_z":5.0,"voxel_h":0.05})");
	EXPECT_EQ(c.nx(), 200);
	EXPECT_EQ(c.ny(), 200);
	EXPECT_EQ(c.nz(), 100);
	EXPECT_EQ(c.cell_count(), 200LL * 200LL * 100LL);
}

TEST(Config, OverridesApplied)
{
	Config c = Config::from_json_string(R"({"name":"v1","U":2.5,"d50":0.001,"nu":1.05e-6})");
	EXPECT_EQ(c.name, "v1");
	EXPECT_DOUBLE_EQ(c.U, 2.5);
	EXPECT_DOUBLE_EQ(c.d50, 0.001);
	EXPECT_DOUBLE_EQ(c.nu, 1.05e-6);
}

TEST(Config, RoundTripJson)
{
	Config a = Config::from_json_string(R"({"name":"rt","U":1.5,"rho":1025.0})");
	Config b = Config::from_json_string(a.to_json_string());
	EXPECT_EQ(b.name, "rt");
	EXPECT_DOUBLE_EQ(b.U, 1.5);
	EXPECT_DOUBLE_EQ(b.rho, 1025.0);
	EXPECT_DOUBLE_EQ(b.nu, a.nu);
}

TEST(Config, MalformedJsonThrows)
{
	EXPECT_THROW(Config::from_json_string("{ this is not json"), std::exception);
}

TEST(Config, UnknownKeyThrows)
{
	EXPECT_THROW(Config::from_json_string(R"({"bogus_key":1.0})"), std::runtime_error);
}

TEST(Config, NonPositiveValueThrows)
{
	EXPECT_THROW(Config::from_json_string(R"({"U":-1.0})"), std::runtime_error);
	EXPECT_THROW(Config::from_json_string(R"({"voxel_h":0.0})"), std::runtime_error);
}

TEST(Config, WrongTypeThrows)
{
	EXPECT_THROW(Config::from_json_string(R"({"U":"fast"})"), std::runtime_error);
}

TEST(Config, PorosityRangeChecked)
{
	EXPECT_THROW(Config::from_json_string(R"({"porosity":1.5})"), std::runtime_error);
}
