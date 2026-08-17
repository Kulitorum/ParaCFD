#include "core/geometry/naca_step.h"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
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
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <filesystem>
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
			double p, double t)
		{
			// -0.1036 gives an exactly closed trailing edge. The more familiar
			// -0.1015 coefficient intentionally leaves a finite trailing-edge gap.
			const double yt = 5.0 * t * (0.2969 * std::sqrt(std::max(0.0, x))
				- 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
				- 0.1036 * x * x * x * x);
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
			double m = 0.0, p = 0.0, t = 0.0;
			if (!parse_code(options.four_digit_code, m, p, t, error)) return false;
			if (!(options.chord > 0.0) || !(options.span > 0.0))
			{
				set_error(error, "NACA chord and span must both be positive");
				return false;
			}
			if (options.points_per_side < 17 || options.points_per_side > 4097)
			{
				set_error(error, "NACA points per side must be in [17, 4097]");
				return false;
			}
			if (options.output_path.empty())
			{
				set_error(error, "NACA STEP output path is empty");
				return false;
			}

			std::vector<SectionPoint> upper, lower;
			upper.reserve(options.points_per_side); lower.reserve(options.points_per_side);
			for (int index = options.points_per_side - 1; index >= 0; --index)
			{
				const double beta = kPi * index / (options.points_per_side - 1);
				const double x = 0.5 * (1.0 - std::cos(beta));
				auto sides = naca_sides(x, m, p, t);
				sides.first.x *= options.chord; sides.first.z *= options.chord;
				upper.push_back(sides.first); // trailing edge -> leading edge
			}
			for (int index = 0; index < options.points_per_side; ++index)
			{
				const double beta = kPi * index / (options.points_per_side - 1);
				const double x = 0.5 * (1.0 - std::cos(beta));
				auto sides = naca_sides(x, m, p, t);
				sides.second.x *= options.chord; sides.second.z *= options.chord;
				lower.push_back(sides.second); // leading edge -> trailing edge
			}

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
}
