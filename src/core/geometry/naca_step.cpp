#include "core/geometry/naca_step.h"
#include "core/geometry/naca_theory.h"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pnt.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		constexpr double kMillimetresPerMetre = 1000.0;
		constexpr double kPi = 3.1415926535897932384626433832795;

		void set_error(std::string* error, const std::string& message)
		{
			if (error) *error = message;
		}

		struct SectionPoint { double x = 0.0, z = 0.0; };

		bool parse_code(const std::string& code, double& maximum_camber,
			double& camber_position, double& thickness, std::string* error)
		{
			if (code.size() != 4 || !std::all_of(code.begin(), code.end(),
				[](unsigned char value) { return std::isdigit(value) != 0; }))
			{
				set_error(error, "NACA code must contain exactly four digits (for example 2412)");
				return false;
			}
			maximum_camber = static_cast<double>(code[0] - '0') / 100.0;
			camber_position = static_cast<double>(code[1] - '0') / 10.0;
			thickness = static_cast<double>(10 * (code[2] - '0') + code[3] - '0') / 100.0;
			if (!(thickness > 0.0))
			{
				set_error(error, "NACA thickness must be greater than zero");
				return false;
			}
			if (maximum_camber > 0.0 && !(camber_position > 0.0 && camber_position < 1.0))
			{
				set_error(error, "a cambered four-digit NACA section needs a non-zero camber position digit");
				return false;
			}
			return true;
		}

		std::pair<SectionPoint, SectionPoint> naca_sides(double x, double m,
			double p, double t, NacaGeometryConvention convention)
		{
			double yt = 0.0;
			if (convention == NacaGeometryConvention::nasa_tmr_naca0012)
			{
				// NASA TMR 2DN00 definition.  These normalized coefficients include
				// the report's extension/rescaling of the original NACA 0012 so its
				// sharp trailing edge is exactly at x/c=1.
				yt = 0.594689181 * (0.298222773 * std::sqrt(std::max(0.0, x))
					- 0.127125232 * x - 0.357907906 * x * x
					+ 0.291984971 * x * x * x - 0.105174606 * x * x * x * x);
			}
			else
			{
				const double fourth = convention == NacaGeometryConvention::closed_sharp
					? -0.1036 : -0.1015;
				yt = 5.0 * t * (0.2969 * std::sqrt(std::max(0.0, x))
					- 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
					+ fourth * x * x * x * x);
			}
			double yc = 0.0, slope = 0.0;
			if (m > 0.0 && x < p)
			{
				yc = m * (2.0 * p * x - x * x) / (p * p);
				slope = 2.0 * m * (p - x) / (p * p);
			}
			else if (m > 0.0)
			{
				const double q = 1.0 - p;
				yc = m * (1.0 - 2.0 * p + 2.0 * p * x - x * x) / (q * q);
				slope = 2.0 * m * (p - x) / (q * q);
			}
			const double angle = std::atan(slope), sine = std::sin(angle), cosine = std::cos(angle);
			return {{x - yt * sine, yc + yt * cosine},
				{x + yt * sine, yc - yt * cosine}};
		}

		bool build_section_points(const std::string& code, NacaGeometryConvention convention, double chord,
			int points_per_side, std::vector<SectionPoint>& upper,
			std::vector<SectionPoint>& lower, std::string* error)
		{
			double m = 0.0, p = 0.0, t = 0.0;
			if (!parse_code(code, m, p, t, error)) return false;
			if (convention == NacaGeometryConvention::nasa_tmr_naca0012 && code != "0012")
			{
				set_error(error, "NASA TMR geometry is defined only for NACA 0012");
				return false;
			}
			if (convention == NacaGeometryConvention::nrel_measured_naca4415 && code != "4415")
			{
				set_error(error, "NREL measured geometry is defined only for NACA 4415");
				return false;
			}
			if (!(chord > 0.0))
			{
				set_error(error, "NACA chord must be positive");
				return false;
			}
			if (points_per_side < 17 || points_per_side > 4097)
			{
				set_error(error, "NACA points per side must be in [17, 4097]");
				return false;
			}

			upper.clear(); lower.clear();
			if (convention == NacaGeometryConvention::nrel_measured_naca4415)
			{
				// Table A1 publishes raw inches against an 18-inch desired chord.
				// Normalize only here so the shared source data remain an exact,
				// auditable transcription. The measured trailing edge is at
				// 18.054 inches, so it intentionally remains at x/c=1.003 when the
				// requested chord represents the experiment's nominal 18 inches.
				const double scale = chord / nrel_naca4415_nominal_chord_inches;
				upper.reserve(nrel_naca4415_upper_coordinates_inches.size());
				lower.reserve(nrel_naca4415_lower_coordinates_inches.size());
				for (auto it = nrel_naca4415_upper_coordinates_inches.rbegin();
					it != nrel_naca4415_upper_coordinates_inches.rend(); ++it)
					upper.push_back({it->chord_station * scale, it->ordinate * scale});
				for (const NacaMeasuredCoordinateInches& point
					: nrel_naca4415_lower_coordinates_inches)
					lower.push_back({point.chord_station * scale, point.ordinate * scale});
				return true;
			}
			upper.reserve(points_per_side); lower.reserve(points_per_side);
			for (int index = points_per_side - 1; index >= 0; --index)
			{
				const double beta = kPi * index / (points_per_side - 1);
				const double x = 0.5 * (1.0 - std::cos(beta));
				auto sides = naca_sides(x, m, p, t, convention);
				sides.first.x *= chord; sides.first.z *= chord;
				upper.push_back(sides.first);
			}
			for (int index = 0; index < points_per_side; ++index)
			{
				const double beta = kPi * index / (points_per_side - 1);
				const double x = 0.5 * (1.0 - std::cos(beta));
				auto sides = naca_sides(x, m, p, t, convention);
				sides.second.x *= chord; sides.second.z *= chord;
				lower.push_back(sides.second);
			}
			return true;
		}

		Handle(Geom_BSplineCurve) interpolate(const std::vector<SectionPoint>& points,
			double y_mm)
		{
			Handle(TColgp_HArray1OfPnt) array = new TColgp_HArray1OfPnt(1,
				static_cast<Standard_Integer>(points.size()));
			for (Standard_Integer index = 1; index <= array->Length(); ++index)
			{
				const SectionPoint& point = points[static_cast<std::size_t>(index - 1)];
				array->SetValue(index, gp_Pnt(point.x * kMillimetresPerMetre,
					y_mm, point.z * kMillimetresPerMetre));
			}
			GeomAPI_Interpolate interpolation(array, false, 1e-7);
			interpolation.Perform();
			return interpolation.IsDone() ? interpolation.Curve() : Handle(Geom_BSplineCurve){};
		}
	}

	bool write_naca_4digit_step(const NacaStepOptions& options, std::string* error)
	{
		try
		{
			if (!(options.chord > 0.0) || !(options.span > 0.0))
			{
				set_error(error, "NACA chord and span must both be positive");
				return false;
			}
			if (options.output_path.empty())
			{
				set_error(error, "NACA STEP output path is empty");
				return false;
			}

			std::vector<SectionPoint> upper, lower;
			if (!build_section_points(options.four_digit_code, options.geometry, options.chord,
				options.points_per_side, upper, lower, error)) return false;

			const double y0 = -0.5 * options.span * kMillimetresPerMetre;
			const Handle(Geom_BSplineCurve) upper_curve = interpolate(upper, y0);
			const Handle(Geom_BSplineCurve) lower_curve = interpolate(lower, y0);
			if (upper_curve.IsNull() || lower_curve.IsNull())
			{
				set_error(error, "OpenCascade could not interpolate the NACA section curves");
				return false;
			}
			BRepBuilderAPI_MakeEdge upper_edge(upper_curve), lower_edge(lower_curve);
			if (!upper_edge.IsDone() || !lower_edge.IsDone())
			{
				set_error(error, "OpenCascade could not build NACA section edges");
				return false;
			}
			BRepBuilderAPI_MakeWire wire_builder;
			wire_builder.Add(upper_edge.Edge()); wire_builder.Add(lower_edge.Edge());
			if (options.geometry == NacaGeometryConvention::original_finite_trailing_edge
				|| options.geometry == NacaGeometryConvention::nrel_measured_naca4415)
			{
				const SectionPoint& lower_te = lower.back();
				const SectionPoint& upper_te = upper.front();
				BRepBuilderAPI_MakeEdge trailing_edge(
					gp_Pnt(lower_te.x * kMillimetresPerMetre, y0, lower_te.z * kMillimetresPerMetre),
					gp_Pnt(upper_te.x * kMillimetresPerMetre, y0, upper_te.z * kMillimetresPerMetre));
				if (!trailing_edge.IsDone())
				{
					set_error(error, "OpenCascade could not close the finite NACA trailing edge");
					return false;
				}
				wire_builder.Add(trailing_edge.Edge());
			}
			if (!wire_builder.IsDone())
			{
				set_error(error, "OpenCascade could not close the NACA section wire");
				return false;
			}
			const TopoDS_Wire wire = wire_builder.Wire();
			BRepBuilderAPI_MakeFace face_builder(wire);
			if (!face_builder.IsDone())
			{
				set_error(error, "OpenCascade could not form a face from the NACA section");
				return false;
			}
			BRepPrimAPI_MakePrism prism(face_builder.Face(),
				gp_Vec(0.0, options.span * kMillimetresPerMetre, 0.0));
			prism.Build();
			if (!prism.IsDone() || prism.Shape().IsNull())
			{
				set_error(error, "OpenCascade could not extrude the NACA section");
				return false;
			}

			const std::filesystem::path output(options.output_path);
			if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
			Interface_Static::SetCVal("write.step.unit", "MM");
			STEPControl_Writer writer;
			if (writer.Transfer(prism.Shape(), STEPControl_AsIs) != IFSelect_RetDone
				|| writer.Write(options.output_path.c_str()) != IFSelect_RetDone)
			{
				set_error(error, "OpenCascade STEP writer failed for '" + options.output_path + "'");
				return false;
			}
			if (error) error->clear();
			return true;
		}
		catch (const Standard_Failure& failure)
		{
			set_error(error, std::string("OpenCascade NACA generation failed: ")
				+ (failure.GetMessageString() ? failure.GetMessageString() : "unknown failure"));
			return false;
		}
		catch (const std::exception& exception)
		{
			set_error(error, std::string("NACA generation failed: ") + exception.what());
			return false;
		}
	}

	bool write_naca_4digit_coordinates(const NacaCoordinateOptions& options,
		std::string* error)
	{
		try
		{
			if (options.output_path.empty())
			{
				set_error(error, "NACA coordinate output path is empty");
				return false;
			}
			std::vector<SectionPoint> upper, lower;
			if (!build_section_points(options.four_digit_code, options.geometry, options.chord,
				options.points_per_side, upper, lower, error)) return false;

			const std::filesystem::path output(options.output_path);
			if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
			std::ofstream stream(output, std::ios::trunc);
			if (!stream)
			{
				set_error(error, "could not open NACA coordinate output '" + options.output_path + "'");
				return false;
			}
			stream << "NACA " << options.four_digit_code << ' '
				<< naca_geometry_convention_name(options.geometry) << '\n';
			stream << std::setprecision(16);
			for (const SectionPoint& point : upper) stream << point.x << ' ' << point.z << '\n';
			// Do not duplicate the nose shared by the upper and lower splines.
			for (std::size_t index = 1; index < lower.size(); ++index)
				stream << lower[index].x << ' ' << lower[index].z << '\n';
			if (!stream)
			{
				set_error(error, "failed while writing NACA coordinates to '" + options.output_path + "'");
				return false;
			}
			if (error) error->clear();
			return true;
		}
		catch (const std::exception& exception)
		{
			set_error(error, std::string("NACA coordinate generation failed: ") + exception.what());
			return false;
		}
	}

	bool write_closed_cylinder_step(const ClosedCylinderStepOptions& options,
		std::string* error)
	{
		try
		{
			if (options.output_path.empty())
			{
				set_error(error, "closed-cylinder STEP output path is empty");
				return false;
			}
			if (!(options.radius > 0.0) || !(options.height > 0.0))
			{
				set_error(error, "closed-cylinder radius and height must be positive");
				return false;
			}

			const gp_Ax2 axis(gp_Pnt(options.centre_x * kMillimetresPerMetre,
				options.centre_y * kMillimetresPerMetre,
				options.base_z * kMillimetresPerMetre), gp_Dir(0.0, 0.0, 1.0));
			BRepPrimAPI_MakeCylinder cylinder(axis,
				options.radius * kMillimetresPerMetre,
				options.height * kMillimetresPerMetre);
			cylinder.Build();
			if (!cylinder.IsDone() || cylinder.Shape().IsNull())
			{
				set_error(error, "OpenCascade could not construct the closed cylinder");
				return false;
			}

			const std::filesystem::path output(options.output_path);
			if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
			Interface_Static::SetCVal("write.step.unit", "MM");
			STEPControl_Writer writer;
			if (writer.Transfer(cylinder.Shape(), STEPControl_AsIs) != IFSelect_RetDone
				|| writer.Write(options.output_path.c_str()) != IFSelect_RetDone)
			{
				set_error(error, "OpenCascade STEP writer failed for '" + options.output_path + "'");
				return false;
			}
			if (error) error->clear();
			return true;
		}
		catch (const Standard_Failure& failure)
		{
			set_error(error, std::string("OpenCascade closed-cylinder generation failed: ")
				+ (failure.GetMessageString() ? failure.GetMessageString() : "unknown failure"));
			return false;
		}
		catch (const std::exception& exception)
		{
			set_error(error, std::string("closed-cylinder generation failed: ") + exception.what());
			return false;
		}
	}

	const char* naca_geometry_convention_name(NacaGeometryConvention convention)
	{
		switch (convention)
		{
		case NacaGeometryConvention::closed_sharp: return "ParaCFD closed sharp trailing edge";
		case NacaGeometryConvention::original_finite_trailing_edge: return "original finite trailing edge";
		case NacaGeometryConvention::nasa_tmr_naca0012: return "NASA TMR exact NACA 0012";
		case NacaGeometryConvention::nrel_measured_naca4415: return "NREL measured NACA 4415";
		}
		return "unknown NACA geometry";
	}
}
