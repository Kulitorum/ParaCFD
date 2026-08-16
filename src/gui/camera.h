// Slicer-style orbit/pan/zoom camera for the paraglider viewer. Self-contained
// (Qt math only), modelled on cobod-slicer's Render::Camera: azimuth/elevation/distance
// about a target, Z-up. Produces a view-projection QMatrix4x4 for the slice shader.
// Units: metres (SI), matching the simulation domain.
#pragma once

#include <QMatrix4x4>
#include <QVector3D>

#include <algorithm>
#include <cmath>

namespace paracfd::gui
{
	class Camera
	{
	public:
		void setTarget(const QVector3D& t) { target_ = t; }
		void setDistance(float d) { distance_ = std::max(1e-3f, d); }
		void setViewport(int w, int h) { width_ = std::max(1, w); height_ = std::max(1, h); }

		// Frame a domain [0,Lx]x[0,Ly]x[0,Lz]: target its centre, back off to fit.
		void frameDomain(float Lx, float Ly, float Lz)
		{
			frameBounds(QVector3D(0.0f, 0.0f, 0.0f), QVector3D(Lx, Ly, Lz), 1.6f);
		}

		// Frame an arbitrary world-space box. This is used for a close wing view while
		// frameDomain remains available for inspecting the asymmetric wake volume.
		void frameBounds(const QVector3D& lo, const QVector3D& hi, float distance_scale = 1.35f)
		{
			target_ = 0.5f * (lo + hi);
			const QVector3D extent = hi - lo;
			const float diag = std::sqrt(QVector3D::dotProduct(extent, extent));
			distance_ = std::max(1e-2f, distance_scale * diag);
			azimuth_ = 35.0f;
			elevation_ = 22.0f;
			panx_ = pany_ = 0.0f;
		}

		// Orbit by mouse deltas (degrees ~ pixels). Elevation clamped off the poles.
		void orbit(float dAz, float dEl)
		{
			azimuth_ += dAz;
			elevation_ = std::clamp(elevation_ + dEl, -89.0f, 89.0f);
		}

		// Pan in the view plane; deltas in pixels, scaled by distance so it tracks.
		void pan(float dxPix, float dyPix)
		{
			float k = distance_ * 0.0015f;
			panx_ += -dxPix * k;
			pany_ += dyPix * k;
		}

		// Zoom: factor > 1 moves closer (wheel up). Exponential, clamped.
		void zoom(float factor)
		{
			distance_ = std::clamp(distance_ / std::max(1e-3f, factor), 1e-2f, 1e5f);
		}

		QVector3D eye() const
		{
			float az = azimuth_ * float(M_PI) / 180.0f;
			float el = elevation_ * float(M_PI) / 180.0f;
			QVector3D dir(std::cos(el) * std::cos(az), std::cos(el) * std::sin(az), std::sin(el));
			return effectiveTarget() + dir * distance_;
		}

		// Unit view direction (eye → target) and the point the camera looks at, both in world
		// metres. Used by the clip plane's "face camera" mode (normal = forward, swept about the
		// target) and its axis-mode auto-side pick (which half is toward the eye).
		QVector3D forward() const
		{
			float az = azimuth_ * float(M_PI) / 180.0f;
			float el = elevation_ * float(M_PI) / 180.0f;
			QVector3D dir(std::cos(el) * std::cos(az), std::cos(el) * std::sin(az), std::sin(el));
			return -dir; // eye is at target + dir*distance, so eye→target is -dir (already unit)
		}
		QVector3D targetPoint() const { return effectiveTarget(); }
		float distance() const { return distance_; }

		QMatrix4x4 view() const
		{
			QMatrix4x4 m;
			m.lookAt(eye(), effectiveTarget(), QVector3D(0, 0, 1));
			return m;
		}

		QMatrix4x4 projection() const
		{
			QMatrix4x4 m;
			float aspect = (float)width_ / (float)height_;
			m.perspective(45.0f, aspect, std::max(1e-3f, distance_ * 0.01f), distance_ * 100.0f);
			return m;
		}

		QMatrix4x4 viewProjection() const { return projection() * view(); }

	private:
		// Pan offset applied in the camera's right/up plane.
		QVector3D effectiveTarget() const
		{
			float az = azimuth_ * float(M_PI) / 180.0f;
			QVector3D right(-std::sin(az), std::cos(az), 0.0f);
			QVector3D up(0, 0, 1);
			return target_ + right * panx_ + up * pany_;
		}

		QVector3D target_{ 0, 0, 0 };
		float distance_ = 10.0f;
		float azimuth_ = 35.0f;   // degrees about +Z
		float elevation_ = 22.0f; // degrees above x-y plane
		float panx_ = 0.0f, pany_ = 0.0f;
		int width_ = 800, height_ = 600;
	};
}
