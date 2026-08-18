// naca_theory.h -- OpenCascade-free analytical and validation references for NACA sections.
#pragma once

#include <array>
#include <cstddef>
#include <string>

namespace paracfd::core
{
struct NacaFourDigitDefinition
{
  std::string code;
  double maximum_camber = 0.0;  // fraction of chord
  double camber_position = 0.0; // fraction of chord
  double thickness = 0.0;       // fraction of chord
};

bool parse_naca_four_digit(const std::string &code, NacaFourDigitDefinition &definition,
    std::string *error = nullptr);
double naca_camber_slope(const NacaFourDigitDefinition &definition, double x_over_chord);
// Classical thin-airfoil zero-lift angle. This is an analytical reference, not an
// empirical viscous polar: finite thickness/Reynolds effects need separate tolerance.
double naca_thin_airfoil_zero_lift_angle(const NacaFourDigitDefinition &definition);
double naca_thin_airfoil_lift_coefficient(const NacaFourDigitDefinition &definition,
    double angle_of_attack_radians);

struct NacaXfoilLiftPoint
{
  double angle_degrees = 0.0;
  double inviscid_cl = 0.0;
  double viscous_cl = 0.0; // free transition, Ncrit from NacaXfoilLiftPolar
  // Boundary-layer transition forced at x/c=0.01 on both surfaces. This is a
  // diagnostic XFOIL prediction, not experimental truth.
  double forced_turbulent_cl = 0.0;
};

struct NacaXfoilLiftPolar
{
  std::string code;
  double reynolds = 666667.0;
  double mach = 0.029;
  double ncrit = 9.0;
  std::array<NacaXfoilLiftPoint, 3> points{};
};

// Independently generated with official MIT XFOIL 6.99. XFOIL loaded the
// exact closed-trailing-edge coordinates emitted by naca_step_generator,
// re-paneled them to 160 panels, then evaluated the attached positive-lift
// regime at 0, 2, and 4 degrees.  Reverse-loaded, highly cambered sections
// are retained as a separate separation diagnostic rather than being used to
// define an accuracy gate. Experimental benchmark data below are authoritative;
// XFOIL is retained only for code-to-code context.
// The table retains three distinct models: inviscid, free-transition viscous,
// and viscous with transition forced at x/c=0.01 on both surfaces.  Callers
// must select the column matching the boundary treatment under validation.
const NacaXfoilLiftPolar *naca_xfoil_lift_polar(const std::string &code);

struct NacaExperimentalLiftPoint
{
  double angle_degrees = 0.0;
  double cl = 0.0;
  double pressure_drag_coefficient = 0.0;
  double moment_coefficient = 0.0;
};

struct NacaExperimentalLiftPolar
{
  const char *id = nullptr;
  const char *profile = nullptr;
  const char *source_url = nullptr;
  const char *geometry = nullptr;
  const char *transition = nullptr;
  double reynolds = 0.0;
  double mach = 0.0;
  const NacaExperimentalLiftPoint *points = nullptr;
  std::size_t point_count = 0;
};

// Raw model-coordinate measurements from Table A1 of the Ohio State/NREL
// NACA 4415 7x10-ft report.  Keeping the published inch values, rather than a
// second rounded normalized copy, makes the CAD generator and geometry gate
// traceable to the same source table.  The two surfaces have different station
// counts and both include the measured leading/trailing-edge points.
struct NacaMeasuredCoordinateInches
{
  double chord_station = 0.0;
  double ordinate = 0.0;
};

inline constexpr double nrel_naca4415_nominal_chord_inches = 18.0;
inline constexpr double nrel_naca4415_measured_trailing_edge_inches = 18.054;
inline constexpr const char *nrel_naca4415_geometry_source_url =
    "https://www.nlr.gov/docs/libraries/wind-images/n4415-7x10.pdf";
inline constexpr std::array<NacaMeasuredCoordinateInches, 49>
    nrel_naca4415_upper_coordinates_inches{{
        {0.0000, 0.0936}, {0.0018, 0.1404}, {0.0036, 0.1530},
        {0.0072, 0.1764}, {0.0108, 0.1980}, {0.0198, 0.2358},
        {0.0522, 0.3222}, {0.0792, 0.3708}, {0.1368, 0.4518},
        {0.2160, 0.5436}, {0.3168, 0.6426}, {0.3924, 0.7074},
        {0.4806, 0.7740}, {0.6156, 0.8658}, {0.7488, 0.9468},
        {0.9036, 1.0315}, {1.0548, 1.1052}, {1.1718, 1.1592},
        {1.3662, 1.2420}, {1.5768, 1.3266}, {1.7424, 1.3878},
        {1.9062, 1.4436}, {2.0394, 1.4868}, {2.4192, 1.5966},
        {2.9718, 1.7262}, {3.3516, 1.8000}, {3.8106, 1.8720},
        {4.2156, 1.9242}, {4.8942, 1.9872}, {5.8392, 2.0304},
        {6.7266, 2.0340}, {7.5186, 2.0106}, {8.5176, 1.9440},
        {9.3708, 1.8630}, {10.3014, 1.7514}, {12.0690, 1.4778},
        {13.8474, 1.1358}, {14.5332, 0.9810}, {14.9994, 0.8712},
        {15.8652, 0.6498}, {16.3368, 0.5220}, {16.8300, 0.3870},
        {17.1738, 0.2916}, {17.3862, 0.2304}, {17.5482, 0.1854},
        {17.6958, 0.1404}, {17.8974, 0.0810}, {17.9478, 0.0648},
        {18.0540, 0.0342}
    }};
inline constexpr std::array<NacaMeasuredCoordinateInches, 51>
    nrel_naca4415_lower_coordinates_inches{{
        {0.0000, 0.0936}, {0.0018, 0.0414}, {0.0036, 0.0306},
        {0.0054, 0.0144}, {0.0108, -0.0054}, {0.0180, -0.0324},
        {0.0468, -0.1062}, {0.0720, -0.1530}, {0.1332, -0.2304},
        {0.2070, -0.2952}, {0.3060, -0.3618}, {0.3798, -0.4032},
        {0.4680, -0.4446}, {0.6012, -0.4968}, {0.7326, -0.5382},
        {0.8874, -0.5778}, {1.0368, -0.6102}, {1.1052, -0.6318},
        {1.3464, -0.6624}, {1.5570, -0.6912}, {1.7226, -0.7092},
        {1.8882, -0.7236}, {2.0196, -0.7326}, {2.3994, -0.7506},
        {2.9520, -0.7578}, {3.3318, -0.7560}, {3.7908, -0.7470},
        {4.1976, -0.7344}, {4.8780, -0.7074}, {5.8230, -0.6624},
        {6.7122, -0.6174}, {7.4934, -0.5742}, {8.5068, -0.5202},
        {9.3618, -0.4716}, {10.2942, -0.4176}, {12.9384, -0.2682},
        {13.8456, -0.2214}, {14.5368, -0.1890}, {15.0012, -0.1674},
        {15.5052, -0.1440}, {15.8688, -0.1278}, {16.3098, -0.1098},
        {16.8372, -0.0846}, {17.3034, -0.0648}, {17.4636, -0.0576},
        {17.6220, -0.0486}, {17.7048, -0.0450}, {17.7840, -0.0414},
        {17.8362, -0.0378}, {17.9064, -0.0342}, {18.0540, -0.0234}
    }};

// NASA TMR 2DN00 / Ladson NASA-TM-4074, 120-grit transition strips at
// x/c=0.05. The exact modified/rescaled geometry must be generated with
// NacaGeometryConvention::nasa_tmr_naca0012. CD is measured total drag and is
// deliberately not an acceptance quantity until ParaCFD implements skin drag.
const NacaExperimentalLiftPolar &nasa_tmr_naca0012_ladson_120_grit();

// Ohio State/NREL 7x10-ft campaign, clean NACA 4415 at Re ~= 1.01 million.
// The reported Cdp is pressure drag (not total drag), which matches ParaCFD's
// explicitly pressure-only aerodynamic load. The model's measured coordinates
// are documented in the cited report and must be used for a strict geometry gate.
const NacaExperimentalLiftPolar &nrel_naca4415_clean_re1m();

// Linear interpolation inside a measured polar. Returns NaN outside its angle
// range instead of silently extrapolating an experiment.
double experimental_lift_coefficient(const NacaExperimentalLiftPolar &polar,
    double angle_degrees);

struct NacaExperimentalLiftSummary
{
  const char *profile = nullptr;
  const char *source_url = nullptr;
  double reynolds = 0.0;
  double lift_curve_slope_per_degree = 0.0;
  double zero_lift_angle_degrees = 0.0;
  double maximum_lift_coefficient = 0.0;
  // Quarter-chord pitching-moment coefficient evaluated at zero lift,
  // denoted C_m0 in NACA Report 460 Table V.
  double quarter_chord_moment_coefficient_at_zero_lift = 0.0;
};

// Corrected infinite-aspect-ratio values from NACA Report 460, Tables II-V.
// The sections deliberately span symmetric/cambered and thin/thick cases.
// Report 460 states that its physical models were faired from the finite
// analytical trailing-edge ordinate to a rounded zero-ordinate edge. Neither
// ParaCFD's -0.1036 sharp closure nor its straight-capped -0.1015 analytical
// section is therefore an exact replica of those historical test articles.
const NacaExperimentalLiftSummary *naca_tr460_lift_summary(const std::string &code);
} // namespace paracfd::core
