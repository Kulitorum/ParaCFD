#include "core/geometry/naca_theory.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

namespace paracfd::core
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;
void set_error(std::string *error, const std::string &message)
{
  if (error) *error = message;
}
} // namespace

bool parse_naca_four_digit(const std::string &code, NacaFourDigitDefinition &definition, std::string *error)
{
  if (code.size() != 4 || !std::all_of(code.begin(), code.end(), [](unsigned char value)
                              { return std::isdigit(value) != 0; }))
  {
	set_error(error, "NACA code must contain exactly four digits");
	return false;
  }
  NacaFourDigitDefinition parsed;
  parsed.code = code;
  parsed.maximum_camber = static_cast<double>(code[0] - '0') / 100.0;
  parsed.camber_position = static_cast<double>(code[1] - '0') / 10.0;
  parsed.thickness = static_cast<double>(10 * (code[2] - '0') + code[3] - '0') / 100.0;
  if (!(parsed.thickness > 0))
  {
	set_error(error, "NACA thickness must be greater than zero");
	return false;
  }
  if (parsed.maximum_camber > 0 && !(parsed.camber_position > 0 && parsed.camber_position < 1))
  {
	set_error(error, "cambered NACA section needs a non-zero camber position");
	return false;
  }
  definition = parsed;
  if (error) error->clear();
  return true;
}

double naca_camber_slope(const NacaFourDigitDefinition &definition, double x)
{
  x = std::clamp(x, 0.0, 1.0);
  const double m = definition.maximum_camber, p = definition.camber_position;
  if (!(m > 0)) return 0;
  if (x < p) return 2 * m * (p - x) / (p * p);
  const double q = 1 - p;
  return 2 * m * (p - x) / (q * q);
}

double naca_thin_airfoil_zero_lift_angle(const NacaFourDigitDefinition &definition)
{
  // alpha_L0 = -(1/pi) integral_0^pi z'_c(theta) (cos(theta)-1) dtheta,
  // x/c=(1-cos(theta))/2. Midpoint quadrature avoids the LE/TE endpoints and
  // is deterministic to substantially better than one microdegree here.
  constexpr int samples = 65536;
  double integral = 0;
  for (int q = 0; q < samples; ++q)
  {
	const double theta = pi * (q + 0.5) / samples, x = 0.5 * (1 - std::cos(theta));
	integral += naca_camber_slope(definition, x) * (std::cos(theta) - 1);
  }
  return -integral / samples;
}

double naca_thin_airfoil_lift_coefficient(const NacaFourDigitDefinition &definition, double alpha)
{
  return 2 * pi * (alpha - naca_thin_airfoil_zero_lift_angle(definition));
}

const NacaXfoilLiftPolar *naca_xfoil_lift_polar(const std::string &code)
{
  // XFOIL 6.99 commands: LOAD exact.dat; PANE; OPER; MACH 0.029;
  // [VISC 666667; ITER 300]; ALFA ... . Free transition used Ncrit=9;
  // the forced-turbulent column additionally used VPAR; XTR 0.01 0.01.
  // Natural-transition NACA2412 did not converge at +4 degrees after 1000
  // iterations, so that contextual (non-oracle) entry remains explicitly NaN.
  const double unavailable = std::numeric_limits<double>::quiet_NaN();
  static const std::array<NacaXfoilLiftPolar, 4> polars = { { { "0012", 666667.0, .029, 9.0, { { { 0.0, 0.0, .0001, 0.0 }, { 2.0, .2415, .1991, .2094 }, { 4.0, .4828, .4496, .4164 } } } },
	  { "2412", 666667.0, .029, 9.0, { { { 0.0, .2593, .2210, .1996 }, { 2.0, .5006, .4685, .4065 }, { 4.0, .7414, unavailable, .6084 } } } },
	  { "4415", 666667.0, .029, 9.0, { { { 0.0, .5346, .4305, .3497 }, { 2.0, .7814, .6756, .5237 }, { 4.0, 1.0273, .9317, .6899 } } } },
	  { "4112", 666667.0, .029, 9.0, { { { 0.0, .4351, .3921, .3531 }, { 2.0, .6775, .5722, .5561 }, { 4.0, .9192, .7648, .7536 } } } } } };
  for (const NacaXfoilLiftPolar &polar : polars)
	if (polar.code == code) return &polar;
  return nullptr;
}

const NacaExperimentalLiftPolar &nasa_tmr_naca0012_ladson_120_grit()
{
  // Direct transcription of CLCD_Ladson_expdata.dat, 120-grit fixed-transition
  // series. Source and geometry are maintained by NASA's Turbulence Modeling
  // Resource; values are not digitized from a plot.
  static constexpr std::array<NacaExperimentalLiftPoint, 18> points = { {
      { -4.01, -0.4466, 0.00843, 0.0 },
      { -2.12, -0.2425, 0.00789, 0.0 },
      { -0.01, -0.0120, 0.00811, 0.0 },
      {  0.01, -0.0122, 0.00804, 0.0 },
      {  2.15,  0.2236, 0.00823, 0.0 },
      {  4.11,  0.4397, 0.00879, 0.0 },
      {  6.01,  0.6487, 0.00842, 0.0 },
      {  8.08,  0.8701, 0.00995, 0.0 },
      { 10.10,  1.0775, 0.01175, 0.0 },
      { 11.23,  1.1849, 0.01248, 0.0 },
      { 12.13,  1.2720, 0.01282, 0.0 },
      { 13.26,  1.3699, 0.01408, 0.0 },
      { 14.30,  1.4571, 0.01628, 0.0 },
      { 15.27,  1.5280, 0.01790, 0.0 },
      { 16.16,  1.5838, 0.02093, 0.0 },
      { 17.24,  1.6347, 0.02519, 0.0 },
      { 18.18,  1.1886, 0.25194, 0.0 },
      { 19.25,  1.1888, 0.28015, 0.0 }
  } };
  static const NacaExperimentalLiftPolar polar{
      "NASA-TMR-2DN00-Ladson-120grit", "0012",
      "https://tmbwg.github.io/turbmodels/NACA0012_validation/CLCD_Ladson_expdata.dat",
      "NASA TMR modified/rescaled sharp-TE NACA 0012",
      "120-grit strips at x/c=0.05, both surfaces", 6.0e6, 0.15,
      points.data(), points.size() };
  return polar;
}

const NacaExperimentalLiftPolar &nrel_naca4415_clean_re1m()
{
  // Corrected values in N4415c100.txt from the report's raw archive.
  static constexpr std::array<NacaExperimentalLiftPoint, 8> points = { {
      { -4.1, 0.018, -0.00884, -0.0968 },
      { -2.1, 0.226, -0.00608, -0.0955 },
      {  0.0, 0.421, -0.00049, -0.0939 },
      {  2.1, 0.629,  0.00495, -0.0932 },
      {  4.1, 0.858,  0.01259, -0.0947 },
      {  6.2, 1.022,  0.01971, -0.0846 },
      {  8.1, 1.171,  0.02570, -0.0746 },
      { 10.2, 1.261,  0.03605, -0.0598 }
  } };
  static const NacaExperimentalLiftPolar polar{
      "NREL-OSU-NACA4415-clean-Re1m", "4415",
      "https://www.nlr.gov/docs/libraries/wind-docs/zip/n4415.zip?sfvrsn=96f27bf7_1",
      "measured NREL/OSU NACA 4415 model ordinates (report table A1)",
      "clean/free transition", 1.01e6, 0.0,
      points.data(), points.size() };
  return polar;
}

double experimental_lift_coefficient(const NacaExperimentalLiftPolar &polar,
    double angle_degrees)
{
  if (!polar.points || polar.point_count == 0 || !std::isfinite(angle_degrees))
    return std::numeric_limits<double>::quiet_NaN();
  if (angle_degrees < polar.points[0].angle_degrees
      || angle_degrees > polar.points[polar.point_count - 1].angle_degrees)
    return std::numeric_limits<double>::quiet_NaN();
  for (std::size_t q = 0; q < polar.point_count; ++q)
  {
    if (angle_degrees == polar.points[q].angle_degrees) return polar.points[q].cl;
    if (q + 1 < polar.point_count && angle_degrees < polar.points[q + 1].angle_degrees)
    {
      const auto &a = polar.points[q], &b = polar.points[q + 1];
      const double t = (angle_degrees - a.angle_degrees)
          / (b.angle_degrees - a.angle_degrees);
      return a.cl + t * (b.cl - a.cl);
    }
  }
  return polar.points[polar.point_count - 1].cl;
}

const NacaExperimentalLiftSummary *naca_tr460_lift_summary(const std::string &code)
{
  static constexpr const char *source =
      "https://ntrs.nasa.gov/api/citations/19930091108/downloads/19930091108.pdf";
  // Reynolds numbers are read from the individual section figures. Lift
  // slope, zero-lift angle, CLmax, and quarter-chord moment at zero lift are
  // direct entries in Tables II, III, IV, and V respectively. Table V calls
  // the last quantity C_m0; Table VI defines C_m(c/4) = C_m0 + n C_L.
  static const std::array<NacaExperimentalLiftSummary, 5> summaries = { {
      { "0012", source, 3.23e6, 0.101,  0.0, 1.53, -0.002 },
      { "0021", source, 3.19e6, 0.094, -0.1, 1.38, -0.001 },
      { "2412", source, 3.25e6, 0.103, -1.7, 1.59, -0.042 },
      { "4412", source, 3.20e6, 0.100, -3.9, 1.61, -0.087 },
      { "6412", source, 3.09e6, 0.102, -5.7, 1.65, -0.129 }
  } };
  for (const auto &summary : summaries) if (code == summary.profile) return &summary;
  return nullptr;
}
} // namespace paracfd::core
