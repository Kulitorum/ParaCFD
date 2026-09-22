#include "gui/paraglider_sim_worker.h"
#include "core/fluid/amr_sampling.h"

#include <QThread>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace
{
	float robust_surface_range(const std::vector<float>& first,
		const std::vector<float>* second = nullptr)
	{
		std::vector<float> magnitude;
		magnitude.reserve(first.size() + (second ? second->size() : 0));
		auto append = [&](const std::vector<float>& values)
		{
			for (float value : values)
				if (std::isfinite(value)) magnitude.push_back(std::abs(value));
		};
		append(first);if(second)append(*second);
		if(magnitude.empty())return 1e-5f;
		// A few sharp-edge/sliver triangles can legitimately have much larger Cp than
		// the rest of the canopy. Let the upper one percent saturate instead of
		// collapsing all useful surface contrast into the green midpoint.
		const std::size_t index=static_cast<std::size_t>(0.99*(magnitude.size()-1));
		std::nth_element(magnitude.begin(),magnitude.begin()+index,magnitude.end());
		return std::max(magnitude[index],1e-5f);
	}

	float robust_volume_focus(const std::vector<float>& values)
	{
		std::vector<float> magnitude;
		magnitude.reserve(values.size());
		for(float value:values)if(std::isfinite(value))magnitude.push_back(std::abs(value));
		if(magnitude.empty())return 1.0f;
		const std::size_t index=static_cast<std::size_t>(0.98*(magnitude.size()-1));
		std::nth_element(magnitude.begin(),magnitude.begin()+index,magnitude.end());
		return std::max(magnitude[index],1e-8f);
	}

	float robust_positive_focus(const std::vector<float>& values)
	{
		std::vector<float> positive;positive.reserve(values.size());
		for(float value:values)if(std::isfinite(value)&&value>0)positive.push_back(value);
		if(positive.empty())return 1e-8f;
		const std::size_t index=static_cast<std::size_t>(0.98*(positive.size()-1));
		std::nth_element(positive.begin(),positive.begin()+index,positive.end());
		return std::max(positive[index],1e-8f);
	}

	QByteArray compress_as_float(const std::vector<double>& values)
	{
		if(values.empty())return {};
		std::vector<float> packed(values.size());
		std::transform(values.begin(),values.end(),packed.begin(),[](double value)
		{
			return std::isfinite(value)?static_cast<float>(value):0.0f;
		});
		return qCompress(reinterpret_cast<const uchar*>(packed.data()),
			static_cast<qsizetype>(packed.size()*sizeof(float)),1);
	}

	QByteArray compress_as_float_blocks(
		std::initializer_list<const std::vector<double>*> blocks)
	{
		std::size_t count=0;for(const auto* block:blocks)count+=block->size();
		if(!count)return {};
		std::vector<float> packed;packed.reserve(count);
		for(const auto* block:blocks)for(double value:*block)
			packed.push_back(std::isfinite(value)?static_cast<float>(value):0.0f);
		return qCompress(reinterpret_cast<const uchar*>(packed.data()),
			static_cast<qsizetype>(packed.size()*sizeof(float)),1);
	}

	QByteArray compress_float_blocks(
		std::initializer_list<const std::vector<float>*> blocks)
	{
		std::size_t count=0;for(const auto* block:blocks)count+=block->size();
		if(!count)return {};
		QByteArray raw;raw.resize(static_cast<qsizetype>(count*sizeof(float)));
		char* destination=raw.data();
		for(const auto* block:blocks)
		{
			const auto bytes=static_cast<qsizetype>(block->size()*sizeof(float));
			if(bytes){std::memcpy(destination,block->data(),static_cast<std::size_t>(bytes));destination+=bytes;}
		}
		return qCompress(raw,1);
	}

	bool decompress_floats(const QByteArray& compressed,std::size_t count,
		std::vector<double>& values)
	{
		values.clear();if(compressed.isEmpty()||!count)return false;
		const QByteArray raw=qUncompress(compressed);
		if(raw.size()!=static_cast<qsizetype>(count*sizeof(float)))return false;
		values.resize(count);const float* source=reinterpret_cast<const float*>(raw.constData());
		std::transform(source,source+count,values.begin(),[](float value)
		{
			return std::isfinite(value)?static_cast<double>(value):0.0;
		});
		return true;
	}

	bool decompress_surface(const QByteArray& compressed,std::size_t triangles,
		std::vector<float>& plus,std::vector<float>& minus,std::vector<float>& delta,
		std::vector<float>& forces)
	{
		plus.clear();minus.clear();delta.clear();forces.clear();
		if(compressed.isEmpty()||!triangles)return false;
		const std::size_t count=6*triangles;const QByteArray raw=qUncompress(compressed);
		if(raw.size()!=static_cast<qsizetype>(count*sizeof(float)))return false;
		const float* source=reinterpret_cast<const float*>(raw.constData());
		plus.assign(source,source+triangles);source+=triangles;
		minus.assign(source,source+triangles);source+=triangles;
		delta.assign(source,source+triangles);source+=triangles;
		forces.assign(source,source+3*triangles);return true;
	}

	bool field_needs_pressure(paracfd::gui::Field field)
	{
		using paracfd::gui::Field;
		return field==Field::PressureDelta||field==Field::PressureCoefficient||
			field==Field::PressureGradient;
	}

	double coarse_velocity_component(const paracfd::core::MacGrid& grid,
		const std::vector<double>& u,const std::vector<double>& v,
		const std::vector<double>& w,int component,int i,int j,int k)
	{
		i=std::clamp(i,0,grid.nx-1);j=std::clamp(j,0,grid.ny-1);k=std::clamp(k,0,grid.nz-1);
		if(component==0)return 0.5*(u[grid.uidx(i,j,k)]+u[grid.uidx(i+1,j,k)]);
		if(component==1)return 0.5*(v[grid.vidx(i,j,k)]+v[grid.vidx(i,j+1,k)]);
		return 0.5*(w[grid.widx(i,j,k)]+w[grid.widx(i,j,k+1)]);
	}

	double coarse_face_sample(const paracfd::core::MacGrid& grid,
		const std::vector<double>& values,int component,double x,double y,double z)
	{
		double coordinate[3]={x/grid.h-(component==0?0.0:0.5),
			y/grid.h-(component==1?0.0:0.5),z/grid.h-(component==2?0.0:0.5)};
		int base[3];double fraction[3];
		for(int axis=0;axis<3;++axis)
		{
			base[axis]=static_cast<int>(std::floor(coordinate[axis]));
			fraction[axis]=coordinate[axis]-base[axis];
		}
		const int count[3]={component==0?grid.nx+1:grid.nx,
			component==1?grid.ny+1:grid.ny,component==2?grid.nz+1:grid.nz};
		double result=0;
		for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx)
		{
			const int i=std::clamp(base[0]+dx,0,count[0]-1);
			const int j=std::clamp(base[1]+dy,0,count[1]-1);
			const int k=std::clamp(base[2]+dz,0,count[2]-1);
			const double weight=(dx?fraction[0]:1-fraction[0])*
				(dy?fraction[1]:1-fraction[1])*(dz?fraction[2]:1-fraction[2]);
			const std::size_t index=component==0?grid.uidx(i,j,k):
				(component==1?grid.vidx(i,j,k):grid.widx(i,j,k));
			result+=weight*values[index];
		}
		return result;
	}

	double coarse_cell_sample(const paracfd::core::MacGrid& grid,
		const std::vector<double>& values,double x,double y,double z)
	{
		double coordinate[3]={x/grid.h-0.5,y/grid.h-0.5,z/grid.h-0.5};
		int base[3];double fraction[3];
		for(int axis=0;axis<3;++axis)
		{
			base[axis]=static_cast<int>(std::floor(coordinate[axis]));
			fraction[axis]=coordinate[axis]-base[axis];
		}
		double result=0;
		for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx)
		{
			const int i=std::clamp(base[0]+dx,0,grid.nx-1);
			const int j=std::clamp(base[1]+dy,0,grid.ny-1);
			const int k=std::clamp(base[2]+dz,0,grid.nz-1);
			const double weight=(dx?fraction[0]:1-fraction[0])*
				(dy?fraction[1]:1-fraction[1])*(dz?fraction[2]:1-fraction[2]);
			result+=weight*values[grid.pidx(i,j,k)];
		}
		return result;
	}

	float coarse_scalar_at_point(const paracfd::core::MacGrid& grid,
		const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& w,
		const std::vector<double>& pressure,bool have_velocity,bool have_pressure,
		paracfd::gui::Field field,double x,double y,double z,double rho,double reference_speed)
	{
		using paracfd::gui::Field;
		if(field_needs_pressure(field)&&!have_pressure)return std::numeric_limits<float>::quiet_NaN();
		if(!field_needs_pressure(field)&&!have_velocity)return std::numeric_limits<float>::quiet_NaN();
		auto velocity=[&](int component,double px,double py,double pz)
		{
			return coarse_face_sample(grid,component==0?u:(component==1?v:w),component,px,py,pz);
		};
		auto derivative=[&](int component,int axis)
		{
			double lo[3]={x,y,z},hi[3]={x,y,z};const double extent[3]={grid.nx*grid.h,grid.ny*grid.h,grid.nz*grid.h};
			lo[axis]=std::max(0.0,lo[axis]-grid.h);hi[axis]=std::min(extent[axis],hi[axis]+grid.h);
			if(!(hi[axis]>lo[axis]))return 0.0;
			return(velocity(component,hi[0],hi[1],hi[2])-velocity(component,lo[0],lo[1],lo[2]))/
				(hi[axis]-lo[axis]);
		};
		if(field==Field::PressureDelta)return static_cast<float>(coarse_cell_sample(grid,pressure,x,y,z));
		if(field==Field::PressureCoefficient)
		{
			const double dynamic_pressure=0.5*rho*reference_speed*reference_speed;
			return dynamic_pressure>0?static_cast<float>(coarse_cell_sample(grid,pressure,x,y,z)/dynamic_pressure):0.0f;
		}
		if(field==Field::PressureGradient)
		{
			double gradient[3]{};const double extent[3]={grid.nx*grid.h,grid.ny*grid.h,grid.nz*grid.h};
			for(int axis=0;axis<3;++axis)
			{
				double lo[3]={x,y,z},hi[3]={x,y,z};lo[axis]=std::max(0.0,lo[axis]-grid.h);
				hi[axis]=std::min(extent[axis],hi[axis]+grid.h);
				if(hi[axis]>lo[axis])gradient[axis]=
					(coarse_cell_sample(grid,pressure,hi[0],hi[1],hi[2])-
					coarse_cell_sample(grid,pressure,lo[0],lo[1],lo[2]))/(hi[axis]-lo[axis]);
			}
			return static_cast<float>(std::sqrt(gradient[0]*gradient[0]+
				gradient[1]*gradient[1]+gradient[2]*gradient[2]));
		}
		const double components[3]={velocity(0,x,y,z),velocity(1,x,y,z),velocity(2,x,y,z)};
		if(field==Field::VelU)return static_cast<float>(components[0]);
		if(field==Field::VelV)return static_cast<float>(components[1]);
		if(field==Field::VelW)return static_cast<float>(components[2]);
		if(field==Field::SpeedMag)return static_cast<float>(std::sqrt(components[0]*components[0]+
			components[1]*components[1]+components[2]*components[2]));
		double gradient[3][3];
		for(int row=0;row<3;++row)for(int column=0;column<3;++column)
			gradient[row][column]=derivative(row,column);
		if(field==Field::VorticityMagnitude)
		{
			const double ox=gradient[2][1]-gradient[1][2],oy=gradient[0][2]-gradient[2][0],
				oz=gradient[1][0]-gradient[0][1];
			return static_cast<float>(std::sqrt(ox*ox+oy*oy+oz*oz));
		}
		double strain2=0,rotation2=0;
		for(int row=0;row<3;++row)for(int column=0;column<3;++column)
		{
			const double strain=0.5*(gradient[row][column]+gradient[column][row]);
			const double rotation=0.5*(gradient[row][column]-gradient[column][row]);
			strain2+=strain*strain;rotation2+=rotation*rotation;
		}
		return static_cast<float>(0.5*(rotation2-strain2));
	}

	float coarse_scalar(const paracfd::core::MacGrid& grid,
		const std::vector<double>& u,const std::vector<double>& v,const std::vector<double>& w,
		const std::vector<double>& pressure,bool have_velocity,bool have_pressure,
		paracfd::gui::Field field,int i,int j,int k,double rho,double reference_speed)
	{
		using paracfd::gui::Field;
		if(field_needs_pressure(field)&&!have_pressure)return std::numeric_limits<float>::quiet_NaN();
		if(!field_needs_pressure(field)&&!have_velocity)return std::numeric_limits<float>::quiet_NaN();
		auto component=[&](int c,int ii,int jj,int kk)
		{
			return coarse_velocity_component(grid,u,v,w,c,ii,jj,kk);
		};
		auto derivative=[&](int component_index,int axis)
		{
			const int coordinate=axis==0?i:axis==1?j:k;
			const int count=axis==0?grid.nx:axis==1?grid.ny:grid.nz;
			const int a=std::clamp(coordinate-1,0,count-1),b=std::clamp(coordinate+1,0,count-1);
			if(a==b)return 0.0;
			int ia=i,ja=j,ka=k,ib=i,jb=j,kb=k;
			if(axis==0){ia=a;ib=b;}else if(axis==1){ja=a;jb=b;}else{ka=a;kb=b;}
			return(component(component_index,ib,jb,kb)-component(component_index,ia,ja,ka))/((b-a)*grid.h);
		};
		if(field==Field::PressureDelta)return static_cast<float>(pressure[grid.pidx(i,j,k)]);
		if(field==Field::PressureCoefficient)
		{
			const double dynamic_pressure=0.5*rho*reference_speed*reference_speed;
			return dynamic_pressure>0?static_cast<float>(pressure[grid.pidx(i,j,k)]/dynamic_pressure):0.0f;
		}
		if(field==Field::PressureGradient)
		{
			auto pd=[&](int axis)
			{
				const int coordinate=axis==0?i:axis==1?j:k;
				const int count=axis==0?grid.nx:axis==1?grid.ny:grid.nz;
				const int a=std::clamp(coordinate-1,0,count-1),b=std::clamp(coordinate+1,0,count-1);
				if(a==b)return 0.0;
				int ia=i,ja=j,ka=k,ib=i,jb=j,kb=k;
				if(axis==0){ia=a;ib=b;}else if(axis==1){ja=a;jb=b;}else{ka=a;kb=b;}
				return(pressure[grid.pidx(ib,jb,kb)]-pressure[grid.pidx(ia,ja,ka)])/((b-a)*grid.h);
			};
			const double x=pd(0),y=pd(1),z=pd(2);return static_cast<float>(std::sqrt(x*x+y*y+z*z));
		}
		const double x=component(0,i,j,k),y=component(1,i,j,k),z=component(2,i,j,k);
		if(field==Field::VelU)return static_cast<float>(x);
		if(field==Field::VelV)return static_cast<float>(y);
		if(field==Field::VelW)return static_cast<float>(z);
		if(field==Field::SpeedMag)return static_cast<float>(std::sqrt(x*x+y*y+z*z));
		double gradient[3][3];
		for(int row=0;row<3;++row)for(int column=0;column<3;++column)
			gradient[row][column]=derivative(row,column);
		if(field==Field::VorticityMagnitude)
		{
			const double ox=gradient[2][1]-gradient[1][2],oy=gradient[0][2]-gradient[2][0],oz=gradient[1][0]-gradient[0][1];
			return static_cast<float>(std::sqrt(ox*ox+oy*oy+oz*oz));
		}
		double strain2=0,rotation2=0;
		for(int row=0;row<3;++row)for(int column=0;column<3;++column)
		{
			const double strain=0.5*(gradient[row][column]+gradient[column][row]);
			const double rotation=0.5*(gradient[row][column]-gradient[column][row]);
			strain2+=strain*strain;rotation2+=rotation*rotation;
		}
		return static_cast<float>(0.5*(rotation2-strain2));
	}
}

namespace paracfd::gui
{
	ParagliderSimWorker::ParagliderSimWorker(
		std::unique_ptr<paracfd::core::ExternalAeroCore> core)
		: core_(std::move(core))
	{
		display_amr_ = std::make_unique<paracfd::core::AmrHostFields>(core_->hierarchy());
		latest_flow_change_=std::numeric_limits<double>::infinity();
	}

	void ParagliderSimWorker::setPlaying(bool playing)
	{
		if(playing)
		{
			endPlayback();
			// Publish the reset before releasing the worker loop. Otherwise the worker can
			// observe playing=true and evaluate one sample against the previous run's
			// convergence history. A queued manual step must not survive into Play either.
			step_requests_.store(0);auto_paused_.store(false);
			pause_reason_.store(SimulationPauseReason::None);settling_reset_.store(true);
			playing_.store(true);
		}
		else
		{
			playing_.store(false);step_requests_.store(0);
			if(!auto_paused_.load())settling_reset_.store(true);
		}
		std::lock_guard lock(snapshot_mutex_);snapshot_.playing=playing;snapshot_.auto_paused=auto_paused_.load();
		if(playing)
		{
			snapshot_.settling_ready=false;snapshot_.mean_force_ready=false;
			snapshot_.bounded_mean_force_ready=false;snapshot_.bounded_mean_force_complete=false;
			snapshot_.pause_reason=SimulationPauseReason::None;snapshot_.settling_score=0;
			snapshot_.settling_force_drift=0;snapshot_.settling_force_rms=0;
			snapshot_.mean_force=snapshot_.mean_pressure_force=snapshot_.mean_viscous_force={};snapshot_.mean_force_drift=0;snapshot_.mean_force_rms=0;
			snapshot_.bounded_mean_force=snapshot_.bounded_mean_pressure_force=snapshot_.bounded_mean_viscous_force={};snapshot_.bounded_mean_force_coverage=0;
			snapshot_.flow_throughs=0;
		}
		++snapshot_.generation;
	}

	void ParagliderSimWorker::stepOnce()
	{
		// Step is a paused-state operation. Enqueuing it while Play is active leaves a
		// latent request that otherwise executes immediately after an auto-pause and
		// mutates the supposedly frozen convergence endpoint.
		if(playing_.load())return;
		auto_paused_.store(false);pause_reason_.store(SimulationPauseReason::None);settling_reset_.store(true);step_requests_.fetch_add(1);
		std::lock_guard lock(snapshot_mutex_);snapshot_.auto_paused=false;snapshot_.pause_reason=SimulationPauseReason::None;
		snapshot_.settling_ready=false;snapshot_.mean_force_ready=false;
		snapshot_.bounded_mean_force_ready=false;snapshot_.bounded_mean_force_complete=false;
		snapshot_.settling_score=0;snapshot_.settling_force_drift=0;snapshot_.settling_force_rms=0;
		snapshot_.mean_force=snapshot_.mean_pressure_force=snapshot_.mean_viscous_force={};snapshot_.mean_force_drift=0;snapshot_.mean_force_rms=0;
		snapshot_.bounded_mean_force=snapshot_.bounded_mean_pressure_force=snapshot_.bounded_mean_viscous_force={};snapshot_.bounded_mean_force_coverage=0;snapshot_.flow_throughs=0;
		++snapshot_.generation;
	}

	void ParagliderSimWorker::configureAutoPause(bool enabled,double sensitivity)
	{
		auto_pause_enabled_.store(enabled);auto_pause_sensitivity_.store(std::clamp(sensitivity,0.0,1.0));settling_reset_.store(true);
		std::lock_guard lock(snapshot_mutex_);snapshot_.auto_pause_enabled=enabled;++snapshot_.generation;
	}

	void ParagliderSimWorker::configureConvergenceExit(double relative_mean_tolerance,
		double maximum_flow_throughs)
	{
		mean_force_tolerance_.store(std::clamp(relative_mean_tolerance,0.001,0.2));
		maximum_flow_throughs_.store(std::clamp(maximum_flow_throughs,1.0,20.0));
		settling_reset_.store(true);
	}

	void ParagliderSimWorker::configureRecording(const RecordingOptions& options)
	{
		// Configuration is normally supplied before the worker thread starts. Keeping the
		// mutation locked also makes a late UI change deterministic: it starts a fresh dataset.
		endPlayback();
		std::lock_guard lock(recording_mutex_);
		recording_options_=options;recorded_frames_.clear();recorded_compressed_bytes_=0;
		++recording_epoch_;next_recording_time_=0;recording_started_.store(false);
	}

	PlaybackDatasetInfo ParagliderSimWorker::playbackDatasetInfo() const
	{
		std::lock_guard lock(recording_mutex_);PlaybackDatasetInfo info;
		info.channels=recording_options_;info.frame_count=recorded_frames_.size();
		info.compressed_bytes=recorded_compressed_bytes_;
		if(!recorded_frames_.empty())
		{
			info.first_time=recorded_frames_.front().physical_time;
			info.last_time=recorded_frames_.back().physical_time;
		}
		info.recording=recording_started_.load()&&!stop_.load();return info;
	}

	std::size_t ParagliderSimWorker::playbackFrameAtTime(double physical_time,int direction) const
	{
		std::lock_guard lock(recording_mutex_);if(recorded_frames_.empty())return 0;
		if(direction<0)
		{
			const auto frame=std::lower_bound(recorded_frames_.begin(),recorded_frames_.end(),physical_time,
				[](const RecordedFrame& candidate,double time){return candidate.physical_time<time;});
			return frame==recorded_frames_.end()?recorded_frames_.size()-1:
				static_cast<std::size_t>(std::distance(recorded_frames_.begin(),frame));
		}
		const auto frame=std::upper_bound(recorded_frames_.begin(),recorded_frames_.end(),physical_time,
			[](double time,const RecordedFrame& candidate){return time<candidate.physical_time;});
		return frame==recorded_frames_.begin()?0:
			static_cast<std::size_t>(std::distance(recorded_frames_.begin(),frame)-1);
	}

	void ParagliderSimWorker::endPlayback()
	{
		std::lock_guard lock(flow_mutex_);playback_active_.store(false);
		playback_u_.clear();playback_v_.clear();playback_w_.clear();playback_p_.clear();
		playback_has_velocity_=playback_has_pressure_=false;++playback_generation_;
	}

	bool ParagliderSimWorker::setPlaybackFrame(std::size_t index,PlaybackFrameSnapshot& snapshot)
	{
		if(playing_.load())return false;
		RecordedFrame frame;std::size_t frame_count=0;
		{
			std::lock_guard lock(recording_mutex_);frame_count=recorded_frames_.size();
			if(index>=frame_count)return false;frame=recorded_frames_[index];
		}
		std::vector<double> u,v,w,p;bool have_velocity=false,have_pressure=false;
		const std::size_t u_count=static_cast<std::size_t>(frame.grid.u_count());
		const std::size_t v_count=static_cast<std::size_t>(frame.grid.v_count());
		const std::size_t w_count=static_cast<std::size_t>(frame.grid.w_count());
		if(!frame.velocity.isEmpty())
		{
			std::vector<double> packed;
			if(decompress_floats(frame.velocity,u_count+v_count+w_count,packed))
			{
				u.assign(packed.begin(),packed.begin()+static_cast<std::ptrdiff_t>(u_count));
				v.assign(packed.begin()+static_cast<std::ptrdiff_t>(u_count),
					packed.begin()+static_cast<std::ptrdiff_t>(u_count+v_count));
				w.assign(packed.begin()+static_cast<std::ptrdiff_t>(u_count+v_count),packed.end());
				have_velocity=true;
			}
		}
		if(!frame.pressure.isEmpty())have_pressure=decompress_floats(
			frame.pressure,static_cast<std::size_t>(frame.grid.p_count()),p);
		{
			std::lock_guard lock(flow_mutex_);playback_grid_=frame.grid;
			playback_u_=std::move(u);playback_v_=std::move(v);playback_w_=std::move(w);playback_p_=std::move(p);
			playback_has_velocity_=have_velocity;playback_has_pressure_=have_pressure;
			++playback_generation_;playback_active_.store(true);
		}
		snapshot={};snapshot.frame_index=index;snapshot.frame_count=frame_count;
		snapshot.steps=frame.steps;snapshot.physical_time=frame.physical_time;
		snapshot.cp_min=frame.cp_min;snapshot.cp_max=frame.cp_max;
		snapshot.side_cp_min=frame.side_cp_min;snapshot.side_cp_max=frame.side_cp_max;
		snapshot.has_surface_pressure=decompress_surface(frame.surface_pressure,frame.surface_triangles,
			snapshot.cp_plus,snapshot.cp_minus,snapshot.delta_cp,
			snapshot.triangle_pressure_force_xyz);
		return true;
	}

	bool ParagliderSimWorker::withFlowField(
		const std::function<void(const FlowField&)>& fn) const
	{
		std::lock_guard lock(flow_mutex_);
		const bool playback=playback_active_.load();
		const auto& u=playback?playback_u_:flow_u_;
		const auto& v=playback?playback_v_:flow_v_;
		const auto& w=playback?playback_w_:flow_w_;
		const auto& p=playback?playback_p_:flow_p_;
		if(playback?!playback_has_velocity_:(!flow_ready_||u.empty()))return false;
		FlowField field;
		field.u=u.data();field.v=v.data();field.w=w.data();field.p=p.empty()?nullptr:p.data();
		field.grid=playback?playback_grid_:flow_grid_;
		field.generation=playback?playback_generation_:flow_generation_;
		fn(field);
		return true;
	}

	bool ParagliderSimWorker::sampleStateLocked(const paracfd::core::Vec3d& world,
		StateSample& sample,bool need_pressure) const
	{
		using namespace paracfd::core;
		sample={};
		if(!display_amr_ready_)return false;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		const BrickLocation location=hierarchy.locate_finest(world);
		if(!location.found())return false;
		const Vec3d velocity=sample_amr_velocity(hierarchy,*display_amr_,world);
		sample.u=velocity.x;sample.v=velocity.y;sample.w=velocity.z;
		if(need_pressure)
		{
			const int pressure_dof=core_->pressure_system().pressure_dof_at_point(
				core_->embedded_boundary(),world);
			sample.p=pressure_dof>=0&&pressure_dof<static_cast<int>(display_pressure_.size())
				?display_pressure_[pressure_dof]:0.0;
		}
		sample.h=hierarchy.levels()[location.level].h;
		sample.valid=std::isfinite(sample.u)&&std::isfinite(sample.v)&&std::isfinite(sample.w)&&std::isfinite(sample.p);
		return sample.valid;
	}

	bool ParagliderSimWorker::sampleDerivedLocked(const paracfd::core::Vec3d& world,Field field,
		FieldProbeSample& sample,bool all_derived) const
	{
		using namespace paracfd::core;
		const bool sample_pressure=all_derived||field==Field::PressureDelta||
			field==Field::PressureCoefficient||field==Field::PressureGradient;
		StateSample centre;
		if(!sampleStateLocked(world,centre,sample_pressure))return false;
		sample={};sample.valid=true;sample.x=world.x;sample.y=world.y;sample.z=world.z;
		sample.u=centre.u;sample.v=centre.v;sample.w=centre.w;sample.pressure_delta=centre.p;
		const auto& freestream=core_->config().freestream;
		const double dynamic_pressure=0.5*freestream.rho*freestream.speed*freestream.speed;
		sample.pressure_coefficient=dynamic_pressure>0?centre.p/dynamic_pressure:0.0;
		const bool need_pressure=all_derived||field==Field::PressureGradient;
		const bool need_velocity=all_derived||field==Field::VorticityMagnitude||field==Field::QCriterion;
		if(!need_pressure&&!need_velocity)return true;

		const Aabb3d& domain=core_->hierarchy().domain();
		double grad_u[3][3]{};double grad_p[3]{};
		for(int axis=0;axis<3;++axis)
		{
			Vec3d lo=world,hi=world;
			lo[axis]=std::max(domain.lo[axis],world[axis]-centre.h);
			hi[axis]=std::min(std::nextafter(domain.hi[axis],domain.lo[axis]),world[axis]+centre.h);
			StateSample a,b;
			if(!sampleStateLocked(lo,a,need_pressure))a=centre;
			if(!sampleStateLocked(hi,b,need_pressure))b=centre;
			const double distance=hi[axis]-lo[axis];
			if(!(distance>0))continue;
			if(need_velocity)
			{
				grad_u[0][axis]=(b.u-a.u)/distance;
				grad_u[1][axis]=(b.v-a.v)/distance;
				grad_u[2][axis]=(b.w-a.w)/distance;
			}
			if(need_pressure)grad_p[axis]=(b.p-a.p)/distance;
		}
		if(need_pressure)sample.pressure_gradient=std::sqrt(grad_p[0]*grad_p[0]+grad_p[1]*grad_p[1]+grad_p[2]*grad_p[2]);
		if(need_velocity)
		{
			const double ox=grad_u[2][1]-grad_u[1][2];
			const double oy=grad_u[0][2]-grad_u[2][0];
			const double oz=grad_u[1][0]-grad_u[0][1];
			sample.vorticity_magnitude=std::sqrt(ox*ox+oy*oy+oz*oz);
			double strain2=0,rotation2=0;
			for(int row=0;row<3;++row)for(int column=0;column<3;++column)
			{
				const double s=0.5*(grad_u[row][column]+grad_u[column][row]);
				const double o=0.5*(grad_u[row][column]-grad_u[column][row]);
				strain2+=s*s;rotation2+=o*o;
			}
			sample.q_criterion=0.5*(rotation2-strain2);
		}
		return true;
	}

	bool ParagliderSimWorker::sampleAmrSlice(
		const SliceParams& params, std::vector<float>& values, FieldRange& range) const
	{
		if(playback_active_.load())return samplePlaybackSlice(params,values,range);
		using namespace paracfd::core;
		std::lock_guard lock(display_amr_mutex_);
		if (!display_amr_ready_ || params.nu <= 0 || params.nv <= 0) return false;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		const Aabb3d& domain = hierarchy.domain();
		const Vec3d extent = domain.hi - domain.lo;
		values.assign(static_cast<std::size_t>(params.nu) * params.nv, 0.0f);
		range = {};
		range.field_min = std::numeric_limits<float>::max();
		range.field_max = std::numeric_limits<float>::lowest();

		auto inside = [](double value, double lo, double hi)
		{
			return std::clamp(value, lo, std::nextafter(hi, lo));
		};
		for (int b = 0; b < params.nv; ++b)
			for (int a = 0; a < params.nu; ++a)
			{
				float x, y, z;
				slice_vertex_world(params, a, b, x, y, z);
				const Vec3d world{
					inside(domain.lo.x + x, domain.lo.x, domain.lo.x + extent.x),
					inside(domain.lo.y + y, domain.lo.y, domain.lo.y + extent.y),
					inside(domain.lo.z + z, domain.lo.z, domain.lo.z + extent.z)};
				FieldProbeSample sampled;
				if(!sampleDerivedLocked(world,params.field,sampled,false))continue;
				const float speed=static_cast<float>(std::sqrt(sampled.u*sampled.u+sampled.v*sampled.v+sampled.w*sampled.w));
				float scalar = speed;
				switch (params.field)
				{
				case Field::VelU: scalar = static_cast<float>(sampled.u); break;
				case Field::VelV: scalar = static_cast<float>(sampled.v); break;
				case Field::VelW: scalar = static_cast<float>(sampled.w); break;
				case Field::PressureDelta: scalar=static_cast<float>(sampled.pressure_delta);break;
				case Field::PressureCoefficient: scalar=static_cast<float>(sampled.pressure_coefficient);break;
				case Field::PressureGradient: scalar=static_cast<float>(sampled.pressure_gradient);break;
				case Field::VorticityMagnitude: scalar=static_cast<float>(sampled.vorticity_magnitude);break;
				case Field::QCriterion: scalar=static_cast<float>(sampled.q_criterion);break;
				case Field::SpeedMag: break;
				}
				if (!std::isfinite(scalar) || !std::isfinite(speed)) continue;
				values[static_cast<std::size_t>(b) * params.nu + a] = scalar;
				range.field_min = std::min(range.field_min, scalar);
				range.field_max = std::max(range.field_max, scalar);
				range.speed_max = std::max(range.speed_max, speed);
			}
		range.valid = range.field_min <= range.field_max;
		return true;
	}

	bool ParagliderSimWorker::samplePoint(double x,double y,double z,FieldProbeSample& sample) const
	{
		if(playback_active_.load())return samplePlaybackPoint(x,y,z,sample);
		using namespace paracfd::core;
		std::lock_guard lock(display_amr_mutex_);
		if(!display_amr_ready_)return false;
		const Aabb3d& domain=core_->hierarchy().domain();const Vec3d extent=domain.hi-domain.lo;
		if(x<0||y<0||z<0||x>extent.x||y>extent.y||z>extent.z)return false;
		const Vec3d world{std::clamp(domain.lo.x+x,domain.lo.x,std::nextafter(domain.hi.x,domain.lo.x)),
			std::clamp(domain.lo.y+y,domain.lo.y,std::nextafter(domain.hi.y,domain.lo.y)),
			std::clamp(domain.lo.z+z,domain.lo.z,std::nextafter(domain.hi.z,domain.lo.z))};
		if(!sampleDerivedLocked(world,Field::QCriterion,sample,true))return false;
		sample.x=x;sample.y=y;sample.z=z;return true;
	}

	bool ParagliderSimWorker::samplePlaybackSlice(const SliceParams& params,
		std::vector<float>& values,FieldRange& range) const
	{
		paracfd::core::MacGrid grid;std::vector<double> u,v,w,p;
		bool have_velocity=false,have_pressure=false;
		{
			std::lock_guard lock(flow_mutex_);if(!playback_active_.load())return false;
			grid=playback_grid_;have_velocity=playback_has_velocity_;have_pressure=playback_has_pressure_;
			if(have_velocity){u=playback_u_;v=playback_v_;w=playback_w_;}if(have_pressure)p=playback_p_;
		}
		if((field_needs_pressure(params.field)&&!have_pressure)||
			(!field_needs_pressure(params.field)&&!have_velocity)||params.nu<=0||params.nv<=0)return false;
		values.assign(static_cast<std::size_t>(params.nu)*params.nv,0.0f);range={};
		range.field_min=std::numeric_limits<float>::max();range.field_max=std::numeric_limits<float>::lowest();
		for(int row=0;row<params.nv;++row)for(int column=0;column<params.nu;++column)
		{
			float x,y,z;slice_vertex_world(params,column,row,x,y,z);
			const float value=coarse_scalar_at_point(grid,u,v,w,p,have_velocity,have_pressure,
				params.field,x,y,z,params.rho,params.reference_speed);
			if(!std::isfinite(value))continue;
			values[static_cast<std::size_t>(row)*params.nu+column]=value;
			range.field_min=std::min(range.field_min,value);range.field_max=std::max(range.field_max,value);
			if(have_velocity)
			{
				const float speed=coarse_scalar_at_point(grid,u,v,w,p,true,have_pressure,
					Field::SpeedMag,x,y,z,params.rho,params.reference_speed);
				if(std::isfinite(speed))range.speed_max=std::max(range.speed_max,speed);
			}
		}
		range.valid=range.field_min<=range.field_max;return range.valid;
	}

	bool ParagliderSimWorker::samplePlaybackPoint(double x,double y,double z,
		FieldProbeSample& sample) const
	{
		paracfd::core::MacGrid grid;std::vector<double> u,v,w,p;
		{
			std::lock_guard lock(flow_mutex_);if(!playback_active_.load()||
				!playback_has_velocity_||!playback_has_pressure_)return false;
			grid=playback_grid_;u=playback_u_;v=playback_v_;w=playback_w_;p=playback_p_;
		}
		if(x<0||y<0||z<0||x>grid.nx*grid.h||y>grid.ny*grid.h||z>grid.nz*grid.h)return false;
		const auto& freestream=core_->config().freestream;sample={};sample.valid=true;
		sample.x=x;sample.y=y;sample.z=z;
		auto scalar=[&](Field field){return coarse_scalar_at_point(grid,u,v,w,p,true,true,field,
			x,y,z,freestream.rho,freestream.speed);};
		sample.u=scalar(Field::VelU);sample.v=scalar(Field::VelV);sample.w=scalar(Field::VelW);
		sample.pressure_delta=scalar(Field::PressureDelta);
		sample.pressure_coefficient=scalar(Field::PressureCoefficient);
		sample.pressure_gradient=scalar(Field::PressureGradient);
		sample.vorticity_magnitude=scalar(Field::VorticityMagnitude);
		sample.q_criterion=scalar(Field::QCriterion);
		return true;
	}

	bool ParagliderSimWorker::scalarVolume(Field field,ScalarVolume& volume) const
	{
		using namespace paracfd::core;
		MacGrid g;std::vector<double> u_faces,v_faces,w_faces,pressure;std::uint64_t generation=0;
		bool have_velocity=false,have_pressure=false;
		{
			// Keep publication latency bounded: copy the immutable snapshot under the mutex and
			// evaluate derivatives after releasing it. Q and vorticity are intentionally CPU-side
			// visualization products and can otherwise hold this lock for an entire volume pass.
			std::lock_guard lock(flow_mutex_);
			const bool playback=playback_active_.load();
			if(playback)
			{
				have_velocity=playback_has_velocity_;have_pressure=playback_has_pressure_;
				generation=playback_generation_;g=playback_grid_;
				if(have_velocity){u_faces=playback_u_;v_faces=playback_v_;w_faces=playback_w_;}
				if(have_pressure)pressure=playback_p_;
			}
			else
			{
				if(!flow_ready_)return false;have_velocity=!flow_u_.empty();have_pressure=!flow_p_.empty();
				generation=flow_generation_;g=flow_grid_;
				u_faces=flow_u_;v_faces=flow_v_;w_faces=flow_w_;pressure=flow_p_;
			}
			if(volume.valid()&&volume.generation==generation&&volume.field==field)return true;
		}
		if((field_needs_pressure(field)&&!have_pressure)||
			(!field_needs_pressure(field)&&!have_velocity))return false;
		volume={};volume.nx=g.nx;volume.ny=g.ny;volume.nz=g.nz;
		volume.h=static_cast<float>(g.h);volume.generation=generation;volume.field=field;
		volume.values.resize(static_cast<std::size_t>(g.p_count()));
		const auto& freestream=core_->config().freestream;
		volume.minimum=std::numeric_limits<float>::max();volume.maximum=std::numeric_limits<float>::lowest();
		for(int k=0;k<g.nz;++k)for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i)
		{
			float value=coarse_scalar(g,u_faces,v_faces,w_faces,pressure,have_velocity,
				have_pressure,field,i,j,k,freestream.rho,freestream.speed);
			if(!std::isfinite(value))value=0.0f;volume.values[static_cast<std::size_t>(g.pidx(i,j,k))]=value;
			volume.minimum=std::min(volume.minimum,value);volume.maximum=std::max(volume.maximum,value);
		}
		volume.focus_abs=robust_volume_focus(volume.values);volume.focus_positive=robust_positive_focus(volume.values);return volume.valid();
	}

	std::uint64_t ParagliderSimWorker::visualizationGeneration() const
	{
		std::lock_guard lock(flow_mutex_);
		if(playback_active_.load())return playback_generation_;
		return flow_ready_?flow_generation_:0;
	}

	void ParagliderSimWorker::publishFlowField()
	{
		using namespace paracfd::core;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		if (hierarchy.levels().empty()) return;
		std::unique_lock display_lock(display_amr_mutex_);
		core_->download_fields(*display_amr_);
		core_->download_pressure(display_pressure_);
		CompositeCellMomentumState display_momentum;
		if(core_->download_cell_momentum_state(display_momentum))
		{
			// The conservative solver advances cell control-volume momentum. Its projected
			// aperture flux can legitimately be very fast on a vanishingly small opening
			// while carrying negligible volume. Sampling those diagnostic face speeds as a
			// Cartesian volume creates bright grid planes and hides the resolved wake. For
			// display/recording only, reconstruct the MAC field from the evolved momentum;
			// the pressure solve and load integration continue to use the exact apertures.
			CompositeAmrFluxes display_fluxes;
			reconstruct_composite_cell_fluxes_cpu(core_->pressure_system(),display_momentum,
				*display_amr_,display_fluxes);
		}

		const Aabb3d& domain = hierarchy.domain();
		const Vec3d extent = domain.hi - domain.lo;
		// Preserve resolved wake structure when the display budget permits it.
		// Always sampling at level zero discarded finer AMR vortices before
		// arrows, tracers, vorticity and recordings could display them.
		double h = hierarchy.finest_cell_size();
		constexpr double maximum_display_cells = 2.0 * 1024 * 1024;
		while (std::ceil(extent.x / h) * std::ceil(extent.y / h) *
			std::ceil(extent.z / h) > maximum_display_cells) h *= 2;
		MacGrid grid;
		grid.nx = std::max(1, static_cast<int>(std::llround(extent.x / h)));
		grid.ny = std::max(1, static_cast<int>(std::llround(extent.y / h)));
		grid.nz = std::max(1, static_cast<int>(std::llround(extent.z / h)));
		grid.h = h;

		// Resample the continuous composite AMR reconstruction. This copy is solely a
		// visualization/recording product; the CFD fields remain staggered, FP32, and GPU resident.
		const std::size_t cells=static_cast<std::size_t>(grid.p_count());
		std::vector<double> uc(cells),vc(cells),wc(cells),pressure(cells);
		for (int k = 0; k < grid.nz; ++k)
			for (int j = 0; j < grid.ny; ++j)
				for (int i = 0; i < grid.nx; ++i)
				{
					const Vec3d world{domain.lo.x + (i + 0.5) * h,
						domain.lo.y + (j + 0.5) * h, domain.lo.z + (k + 0.5) * h};
					const BrickLocation location=hierarchy.locate_finest(world);
					if(!location.found())continue;
					const std::size_t out=static_cast<std::size_t>(grid.pidx(i,j,k));
					const Vec3d velocity=sample_amr_velocity(hierarchy,*display_amr_,world);
					uc[out]=velocity.x;vc[out]=velocity.y;wc[out]=velocity.z;
					const AmrHostLevelFields& level=display_amr_->levels()[location.level];
					const BrickFieldLayout& layout=level.layout;
					const int bi=location.brick,li=location.cell.x,lj=location.cell.y,lk=location.cell.z;
					const int pressure_dof=core_->pressure_system().pressure_dof_at_point(core_->embedded_boundary(),world);
					pressure[out]=pressure_dof>=0&&pressure_dof<static_cast<int>(display_pressure_.size())
						?display_pressure_[pressure_dof]:static_cast<double>(level.p[layout.cell_index(bi,li,lj,lk)]);
				}
		std::vector<double> u(static_cast<std::size_t>(grid.u_count()));
		std::vector<double> v(static_cast<std::size_t>(grid.v_count()));
		std::vector<double> w(static_cast<std::size_t>(grid.w_count()));
		auto finite_or_zero = [](double value) { return std::isfinite(value) ? value : 0.0; };
		for (int k = 0; k < grid.nz; ++k) for (int j = 0; j < grid.ny; ++j) for (int i = 0; i <= grid.nx; ++i)
		{
			const int left=std::max(0,i-1),right=std::min(grid.nx-1,i);
			u[grid.uidx(i,j,k)]=finite_or_zero(0.5*(uc[grid.pidx(left,j,k)]+uc[grid.pidx(right,j,k)]));
		}
		for (int k = 0; k < grid.nz; ++k) for (int j = 0; j <= grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
		{
			const int below=std::max(0,j-1),above=std::min(grid.ny-1,j);
			v[grid.vidx(i,j,k)]=finite_or_zero(0.5*(vc[grid.pidx(i,below,k)]+vc[grid.pidx(i,above,k)]));
		}
		for (int k = 0; k <= grid.nz; ++k) for (int j = 0; j < grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
		{
			const int back=std::max(0,k-1),front=std::min(grid.nz-1,k);
			w[grid.widx(i,j,k)]=finite_or_zero(0.5*(wc[grid.pidx(i,j,back)]+wc[grid.pidx(i,j,front)]));
		}
		for (double& value : pressure) value = finite_or_zero(value);
		display_amr_ready_ = true;
		display_lock.unlock();

		std::lock_guard lock(flow_mutex_);
		if(flow_u_.size()==u.size()&&flow_v_.size()==v.size()&&flow_w_.size()==w.size())
		{
			double difference2=0,scale2=0;auto accumulate=[&](const std::vector<double>& previous,const std::vector<double>& current){for(std::size_t q=0;q<current.size();++q){const double d=current[q]-previous[q];difference2+=d*d;scale2+=current[q]*current[q];}};accumulate(flow_u_,u);accumulate(flow_v_,v);accumulate(flow_w_,w);latest_flow_change_=std::sqrt(difference2/std::max(scale2,1e-30));
		}
		else latest_flow_change_=std::numeric_limits<double>::infinity();
		flow_grid_ = grid;
		flow_u_.swap(u); flow_v_.swap(v); flow_w_.swap(w); flow_p_.swap(pressure);
		++flow_generation_;
		flow_ready_ = true;
	}

	void ParagliderSimWorker::captureRecordingFrame(double scheduled_time,long long steps)
	{
		RecordingOptions options;std::uint64_t recording_epoch=0;
		{
			std::lock_guard lock(recording_mutex_);options=recording_options_;
			recording_epoch=recording_epoch_;
		}
		if(!options.any())return;
		RecordedFrame frame;frame.physical_time=scheduled_time;frame.steps=steps;
		std::vector<double> u,v,w,p;
		{
			std::lock_guard lock(flow_mutex_);if(!flow_ready_)return;frame.grid=flow_grid_;
			if(options.velocity){u=flow_u_;v=flow_v_;w=flow_w_;}
			if(options.pressure)p=flow_p_;
		}
		if(options.velocity)frame.velocity=compress_as_float_blocks({&u,&v,&w});
		if(options.pressure)frame.pressure=compress_as_float(p);
		if(options.surface_pressure)
		{
			std::vector<float> plus,minus,delta,forces;
			{
				std::lock_guard lock(snapshot_mutex_);plus=snapshot_.cp_plus;minus=snapshot_.cp_minus;
				delta=snapshot_.delta_cp;forces=snapshot_.triangle_pressure_force_xyz;
				frame.cp_min=snapshot_.cp_min;frame.cp_max=snapshot_.cp_max;
				frame.side_cp_min=snapshot_.side_cp_min;frame.side_cp_max=snapshot_.side_cp_max;
			}
			if(!plus.empty()&&plus.size()==minus.size()&&plus.size()==delta.size()&&
				forces.size()==3*plus.size())
			{
				frame.surface_triangles=plus.size();
				frame.surface_pressure=compress_float_blocks({&plus,&minus,&delta,&forces});
			}
		}
		frame.compressed_bytes=static_cast<std::size_t>(frame.velocity.size()+frame.pressure.size()+
			frame.surface_pressure.size());
		{
			std::lock_guard lock(recording_mutex_);
			// A configuration change clears the timeline. Do not append a frame encoded with
			// the superseded channel layout after that clear.
			if(recording_epoch!=recording_epoch_||options.velocity!=recording_options_.velocity||
				options.pressure!=recording_options_.pressure||
				options.surface_pressure!=recording_options_.surface_pressure||
				options.every_step!=recording_options_.every_step)return;
			recorded_compressed_bytes_+=frame.compressed_bytes;recorded_frames_.push_back(std::move(frame));
		}
	}

	void ParagliderSimWorker::captureScheduledRecordingFrames(double physical_time,long long steps)
	{
		for(;;)
		{
			double scheduled=0;
			{
				std::lock_guard lock(recording_mutex_);if(!recording_options_.any()||recording_options_.every_step||
					next_recording_time_>physical_time+1e-10)return;
				scheduled=next_recording_time_;next_recording_time_+=kRecordingInterval;
			}
			captureRecordingFrame(scheduled,steps);
		}
	}

	void ParagliderSimWorker::updateSettling(double physical_time,
		const paracfd::core::AerodynamicLoads& loads,
		const paracfd::core::ExternalAeroConservationStats& conservation)
	{
		const paracfd::core::Vec3d force=loads.viscous_loads_valid?
			loads.total_force:loads.pressure_force;
		if(settling_reset_.exchange(false))
		{
			settling_history_.clear();mean_force_history_.clear();
			mean_pressure_force_history_.clear();mean_viscous_force_history_.clear();
			settling_consecutive_=mean_consecutive_=0;
			settling_score_=settling_force_drift_=settling_force_rms_=flow_throughs_=0;
			settling_epoch_time_=physical_time;settling_ready_=false;
			mean_convergence_={};mean_pressure_convergence_={};mean_viscous_convergence_={};
			bounded_mean_force_=bounded_mean_pressure_force_=bounded_mean_viscous_force_={};
			bounded_mean_force_coverage_=0;
			bounded_mean_force_ready_=bounded_mean_force_complete_=false;
			pause_reason_.store(SimulationPauseReason::None);
		}
		const bool exit_enabled=auto_pause_enabled_.load();
		if(!exit_enabled)return;
		const double speed=std::max(1e-9,std::abs(core_->config().freestream.speed));
		const double flow_time=(core_->hierarchy().domain().hi.x-core_->hierarchy().domain().lo.x)/speed;
		if(!(flow_time>0))return;
		flow_throughs_=std::max(0.0,physical_time-settling_epoch_time_)/flow_time;
		mean_force_history_.push_back({physical_time,force});
		mean_pressure_force_history_.push_back({physical_time,loads.pressure_force});
		mean_viscous_force_history_.push_back({physical_time,
			loads.viscous_loads_valid?loads.viscous_force:paracfd::core::Vec3d{}});
		const double mean_window=std::max(0.5,flow_time),mean_oldest=physical_time-mean_window;
		while(mean_force_history_.size()>2&&mean_force_history_[1].time<mean_oldest)mean_force_history_.pop_front();
		while(mean_pressure_force_history_.size()>2&&mean_pressure_force_history_[1].time<mean_oldest)
			mean_pressure_force_history_.pop_front();
		while(mean_viscous_force_history_.size()>2&&mean_viscous_force_history_[1].time<mean_oldest)
			mean_viscous_force_history_.pop_front();
		mean_convergence_=paracfd::core::assess_aerodynamic_mean_convergence(mean_force_history_,physical_time,mean_window);
		mean_pressure_convergence_=paracfd::core::assess_aerodynamic_mean_convergence(
			mean_pressure_force_history_,physical_time,mean_window);
		mean_viscous_convergence_=paracfd::core::assess_aerodynamic_mean_convergence(
			mean_viscous_force_history_,physical_time,mean_window);
		// Preserve a time-weighted final-window load even if there are not yet enough
		// samples to certify two independent half windows. Maximum-flow exits must not
		// silently turn the sweep's <D>/<L> columns into instantaneous samples.
		const auto bounded_mean=paracfd::core::bounded_aerodynamic_mean(
			mean_force_history_,physical_time,mean_window);
		bounded_mean_force_ready_=bounded_mean.ready;
		bounded_mean_force_complete_=bounded_mean.complete;
		bounded_mean_force_coverage_=bounded_mean.coverage;
		bounded_mean_force_=bounded_mean.force;
		const auto bounded_pressure=paracfd::core::bounded_aerodynamic_mean(
			mean_pressure_force_history_,physical_time,mean_window);
		const auto bounded_viscous=paracfd::core::bounded_aerodynamic_mean(
			mean_viscous_force_history_,physical_time,mean_window);
		bounded_mean_pressure_force_=bounded_pressure.force;
		bounded_mean_viscous_force_=bounded_viscous.force;
		const bool conservative=conservation.volume_weighted_rms_divergence<1e-3;
		if(flow_throughs_>=1.0&&mean_convergence_.ready&&conservative&&
			mean_convergence_.relative_drift<mean_force_tolerance_.load())++mean_consecutive_;
		else mean_consecutive_=0;

		auto finish_run=[&](SimulationPauseReason reason)
		{
			playing_.store(false);step_requests_.store(0);auto_paused_.store(true);
			pause_reason_.store(reason);
			const char* label=reason==SimulationPauseReason::Steady?"steady":
				reason==SimulationPauseReason::MeanConverged?"mean-converged":
				"maximum-flow-throughs";
			std::fprintf(stderr,"[paraglider] auto-pause: reason=%s, t=%.6g s, "
				"flow-throughs=%.3f, settle-score=%.4g, mean-drift=%.4g, mean-rms=%.4g\n",
				label,physical_time,flow_throughs_,settling_score_,
				mean_convergence_.relative_drift,mean_convergence_.current_rms_fraction);
		};
		auto finish_mean_or_maximum=[&]()
		{
			if(mean_consecutive_>=3){finish_run(SimulationPauseReason::MeanConverged);return true;}
			if(flow_throughs_>=maximum_flow_throughs_.load())
			{
				finish_run(SimulationPauseReason::MaximumFlowThroughs);return true;
			}
			return false;
		};
		// A bounded run must still stop when visualization sampling has not yet
		// produced a finite field-change metric. Maximum flow time is an explicit
		// "not converged" result, not a numerical stability guard.
		if(!std::isfinite(latest_flow_change_))
		{
			finish_mean_or_maximum();
			return;
		}
		settling_history_.push_back({physical_time,force,latest_flow_change_});
		const double window=std::max(0.25,0.5*flow_time),oldest=physical_time-window;
		while(!settling_history_.empty()&&settling_history_.front().time<oldest)settling_history_.pop_front();
		if(settling_history_.size()<10)
		{
			finish_mean_or_maximum();
			return;
		}
		const double history_span=physical_time-settling_history_.front().time;
		if(!(history_span>0.1*window)){finish_mean_or_maximum();return;}
		const double analysis_window=std::min(window,history_span);
		const double split=physical_time-0.5*analysis_window;paracfd::core::Vec3d old_mean{},new_mean{};int old_count=0,new_count=0;
		for(const SettlingSample& sample:settling_history_){paracfd::core::Vec3d& mean=sample.time<split?old_mean:new_mean;mean.x+=sample.force.x;mean.y+=sample.force.y;mean.z+=sample.force.z;if(sample.time<split)++old_count;else ++new_count;}
		if(old_count<4||new_count<4){finish_mean_or_maximum();return;}old_mean.x/=old_count;old_mean.y/=old_count;old_mean.z/=old_count;new_mean.x/=new_count;new_mean.y/=new_count;new_mean.z/=new_count;
		auto magnitude=[](const paracfd::core::Vec3d& value){return std::sqrt(value.x*value.x+value.y*value.y+value.z*value.z);};
		const double scale=std::max(1.0,0.5*(magnitude(old_mean)+magnitude(new_mean)));settling_force_drift_=magnitude(new_mean-old_mean)/scale;
		double square=0,flow_change=0;for(const SettlingSample& sample:settling_history_)if(sample.time>=split){const paracfd::core::Vec3d delta=sample.force-new_mean;square+=delta.x*delta.x+delta.y*delta.y+delta.z*delta.z;flow_change+=sample.flow_change;}
		settling_force_rms_=std::sqrt(square/new_count)/scale;flow_change/=new_count;
		const double sensitivity=auto_pause_sensitivity_.load(),range=std::pow(20.0,sensitivity);
		const double drift_tolerance=0.001*range,noise_tolerance=0.003*range,flow_tolerance=0.0005*range;
		settling_score_=std::max({settling_force_drift_/drift_tolerance,settling_force_rms_/noise_tolerance,flow_change/flow_tolerance});settling_ready_=true;
		const bool observation_complete=flow_throughs_>=kAutoPauseMinimumFlowThroughs&&history_span>=0.9*window;
		if(observation_complete&&settling_score_<1.0&&conservative)++settling_consecutive_;else settling_consecutive_=0;

		// Every aerodynamic run retains at least one complete force window, even if
		// the instantaneous-field criterion happens to trigger earlier. Single runs
		// and sweep cases therefore have identical bounded stopping semantics.
		const SimulationPauseReason reason=paracfd::core::select_aerodynamic_run_exit({
			settling_consecutive_,mean_consecutive_,mean_convergence_.ready,
			mean_convergence_.relative_drift,mean_force_tolerance_.load(),
			flow_throughs_,maximum_flow_throughs_.load()});
		if(reason!=SimulationPauseReason::None)finish_run(reason);
	}

	bool ParagliderSimWorker::latestSnapshot(std::uint64_t& generation,
		std::uint64_t& surface_generation, ParagliderDisplaySnapshot& out) const
	{
		std::lock_guard lock(snapshot_mutex_);
		if (snapshot_.generation == generation) return false;
		out.generation = snapshot_.generation;
		out.surface_generation = snapshot_.surface_generation;
		out.steps = snapshot_.steps;
		out.physical_time = snapshot_.physical_time;
		out.dt = snapshot_.dt;
		out.step_ms = snapshot_.step_ms;
		out.projection_ms = snapshot_.projection_ms;
		out.max_diffusion_rate = snapshot_.max_diffusion_rate;
		out.diffusion_substeps = snapshot_.diffusion_substeps;
		out.residual = snapshot_.residual;
		out.max_abs_regular_velocity = snapshot_.max_abs_regular_velocity;
		out.max_abs_special_velocity = snapshot_.max_abs_special_velocity;
		out.max_embedded_cfl_rate = snapshot_.max_embedded_cfl_rate;
		out.effective_cfl = snapshot_.effective_cfl;
		out.pressure_iterations = snapshot_.pressure_iterations;
		out.initialized = snapshot_.initialized;
		out.converged = snapshot_.converged;
		out.conservative_cell_momentum = snapshot_.conservative_cell_momentum;
		out.conservation_valid = snapshot_.conservation_valid;
		out.max_abs_divergence = snapshot_.max_abs_divergence;
		out.volume_weighted_rms_divergence = snapshot_.volume_weighted_rms_divergence;
		out.max_integrated_flux_error = snapshot_.max_integrated_flux_error;
		out.absolute_integrated_flux_error = snapshot_.absolute_integrated_flux_error;
		out.net_integrated_flux_error = snapshot_.net_integrated_flux_error;
		out.flow_change=snapshot_.flow_change;out.settling_score=snapshot_.settling_score;out.settling_force_drift=snapshot_.settling_force_drift;out.settling_force_rms=snapshot_.settling_force_rms;out.flow_throughs=snapshot_.flow_throughs;
		out.mean_force=snapshot_.mean_force;out.mean_pressure_force=snapshot_.mean_pressure_force;out.mean_viscous_force=snapshot_.mean_viscous_force;out.mean_force_drift=snapshot_.mean_force_drift;out.mean_force_rms=snapshot_.mean_force_rms;out.mean_force_ready=snapshot_.mean_force_ready;
		out.bounded_mean_force=snapshot_.bounded_mean_force;out.bounded_mean_pressure_force=snapshot_.bounded_mean_pressure_force;out.bounded_mean_viscous_force=snapshot_.bounded_mean_viscous_force;out.bounded_mean_force_coverage=snapshot_.bounded_mean_force_coverage;out.bounded_mean_force_ready=snapshot_.bounded_mean_force_ready;out.bounded_mean_force_complete=snapshot_.bounded_mean_force_complete;out.pause_reason=snapshot_.pause_reason;
		out.gpu_bytes = snapshot_.gpu_bytes;
		out.playing=snapshot_.playing;out.auto_pause_enabled=snapshot_.auto_pause_enabled;out.auto_paused=snapshot_.auto_paused;out.settling_ready=snapshot_.settling_ready;
		out.error = snapshot_.error;
		out.cp_min = snapshot_.cp_min;
		out.cp_max = snapshot_.cp_max;
		out.pressure_force = snapshot_.pressure_force;
		out.viscous_force = snapshot_.viscous_force;
		out.total_force = snapshot_.total_force;
		out.viscous_loads_valid = snapshot_.viscous_loads_valid;
		out.coefficients_valid = snapshot_.coefficients_valid;
		out.cd_pressure = snapshot_.cd_pressure;
		out.cs_pressure = snapshot_.cs_pressure;
		out.cl_pressure = snapshot_.cl_pressure;
		out.cd_viscous = snapshot_.cd_viscous;
		out.cs_viscous = snapshot_.cs_viscous;
		out.cl_viscous = snapshot_.cl_viscous;
		out.cd = snapshot_.cd;
		out.cs = snapshot_.cs;
		out.cl = snapshot_.cl;
		if (snapshot_.surface_generation != surface_generation)
		{
			out.cp_plus = snapshot_.cp_plus;
			out.cp_minus = snapshot_.cp_minus;
			out.delta_cp = snapshot_.delta_cp;
			out.triangle_pressure_force_xyz=snapshot_.triangle_pressure_force_xyz;
			out.side_cp_min = snapshot_.side_cp_min;
			out.side_cp_max = snapshot_.side_cp_max;
			surface_generation = snapshot_.surface_generation;
		}
		generation = snapshot_.generation;
		return true;
	}

	void ParagliderSimWorker::publish(const paracfd::core::ExternalAeroStepStats& stats,
		bool include_surface,bool update_convergence)
	{
		std::vector<float> cp_plus,cp_minus,delta_cp,triangle_pressure_force_xyz;
		float cp_min = 0.0f, cp_max = 0.0f,side_cp_min=0.0f,side_cp_max=0.0f;
		paracfd::core::AerodynamicLoads loads;
		paracfd::core::ExternalAeroConservationStats conservation;
		const bool have_surface = include_surface && stats.pressure.converged;
		if (have_surface)
		{
			loads = core_->aerodynamic_loads();
			conservation = core_->conservation_stats();
			const float missing=std::numeric_limits<float>::quiet_NaN();
			cp_plus.assign(loads.triangles.size(),missing);
			cp_minus.assign(loads.triangles.size(),missing);
			delta_cp.assign(loads.triangles.size(),missing);
			triangle_pressure_force_xyz.assign(loads.triangles.size()*3,0.0f);
			for (std::size_t triangle = 0; triangle < loads.triangles.size(); ++triangle)
			{
				const auto& source=loads.triangles[triangle];
				if(std::isfinite(source.cp_plus))cp_plus[triangle]=static_cast<float>(source.cp_plus);
				if(std::isfinite(source.cp_minus))cp_minus[triangle]=static_cast<float>(source.cp_minus);
				if(std::isfinite(source.delta_cp))delta_cp[triangle]=static_cast<float>(source.delta_cp);
				triangle_pressure_force_xyz[3*triangle]=static_cast<float>(source.pressure_force.x);
				triangle_pressure_force_xyz[3*triangle+1]=static_cast<float>(source.pressure_force.y);
				triangle_pressure_force_xyz[3*triangle+2]=static_cast<float>(source.pressure_force.z);
			}
			const float maximum=robust_surface_range(delta_cp);
			const float side_maximum=robust_surface_range(cp_plus,&cp_minus);
			cp_min = -maximum;
			cp_max = maximum;
			side_cp_min=-side_maximum;
			side_cp_max=side_maximum;
		}
		if(have_surface&&update_convergence)updateSettling(stats.physical_time,loads,conservation);

		std::lock_guard lock(snapshot_mutex_);
		snapshot_.steps = steps_;
		snapshot_.physical_time = stats.physical_time;
		snapshot_.dt = stats.dt;
		snapshot_.step_ms = stats.gpu_step_ms;
		snapshot_.projection_ms = stats.projection_ms;
		snapshot_.max_diffusion_rate = stats.max_diffusion_rate;
		snapshot_.diffusion_substeps = stats.diffusion_substeps;
		snapshot_.pressure_iterations = stats.pressure.iterations;
		snapshot_.residual = stats.pressure.relative_residual;
		snapshot_.max_abs_regular_velocity = stats.max_abs_regular_velocity;
		snapshot_.max_abs_special_velocity = stats.max_abs_special_velocity;
		snapshot_.max_embedded_cfl_rate = stats.max_embedded_cfl_rate;
		snapshot_.effective_cfl = stats.effective_cfl;
		snapshot_.gpu_bytes = core_->gpu_bytes();
		snapshot_.initialized = core_->initialized();
		snapshot_.converged = stats.pressure.converged;
		snapshot_.conservative_cell_momentum = core_->uses_conservative_cell_momentum();
		snapshot_.playing=playing_.load();snapshot_.auto_pause_enabled=auto_pause_enabled_.load();snapshot_.auto_paused=auto_paused_.load();
		const bool convergence_reset_pending=settling_reset_.load();
		if(convergence_reset_pending)
		{
			snapshot_.settling_ready=false;snapshot_.flow_change=latest_flow_change_;snapshot_.settling_score=0;snapshot_.settling_force_drift=0;snapshot_.settling_force_rms=0;snapshot_.flow_throughs=0;
			snapshot_.mean_force=snapshot_.mean_pressure_force=snapshot_.mean_viscous_force={};snapshot_.mean_force_drift=0;snapshot_.mean_force_rms=0;snapshot_.mean_force_ready=false;
			snapshot_.bounded_mean_force=snapshot_.bounded_mean_pressure_force=snapshot_.bounded_mean_viscous_force={};snapshot_.bounded_mean_force_coverage=0;snapshot_.bounded_mean_force_ready=false;snapshot_.bounded_mean_force_complete=false;snapshot_.pause_reason=SimulationPauseReason::None;
		}
		else
		{
			snapshot_.settling_ready=settling_ready_;snapshot_.flow_change=latest_flow_change_;snapshot_.settling_score=settling_score_;snapshot_.settling_force_drift=settling_force_drift_;snapshot_.settling_force_rms=settling_force_rms_;snapshot_.flow_throughs=flow_throughs_;
			snapshot_.mean_force=(mean_convergence_.previous_mean+mean_convergence_.current_mean)*0.5;snapshot_.mean_force_drift=mean_convergence_.relative_drift;snapshot_.mean_force_rms=mean_convergence_.current_rms_fraction;snapshot_.mean_force_ready=mean_convergence_.ready;
			snapshot_.mean_pressure_force=(mean_pressure_convergence_.previous_mean+mean_pressure_convergence_.current_mean)*0.5;
			snapshot_.mean_viscous_force=(mean_viscous_convergence_.previous_mean+mean_viscous_convergence_.current_mean)*0.5;
			snapshot_.bounded_mean_force=bounded_mean_force_;snapshot_.bounded_mean_pressure_force=bounded_mean_pressure_force_;snapshot_.bounded_mean_viscous_force=bounded_mean_viscous_force_;snapshot_.bounded_mean_force_coverage=bounded_mean_force_coverage_;snapshot_.bounded_mean_force_ready=bounded_mean_force_ready_;snapshot_.bounded_mean_force_complete=bounded_mean_force_complete_;snapshot_.pause_reason=pause_reason_.load();
		}
		snapshot_.error.clear();
		if (have_surface)
		{
			snapshot_.cp_plus = std::move(cp_plus);
			snapshot_.cp_minus = std::move(cp_minus);
			snapshot_.delta_cp = std::move(delta_cp);
			snapshot_.triangle_pressure_force_xyz=std::move(triangle_pressure_force_xyz);
			snapshot_.cp_min = cp_min;
			snapshot_.cp_max = cp_max;
			snapshot_.side_cp_min=side_cp_min;
			snapshot_.side_cp_max=side_cp_max;
			snapshot_.pressure_force = loads.pressure_force;
			snapshot_.viscous_force = loads.viscous_force;
			snapshot_.total_force = loads.total_force;
			snapshot_.viscous_loads_valid = loads.viscous_loads_valid;
			snapshot_.coefficients_valid = loads.force_coefficients_valid;
			snapshot_.cd_pressure = loads.cd_pressure;
			snapshot_.cs_pressure = loads.cs_pressure;
			snapshot_.cl_pressure = loads.cl_pressure;
			snapshot_.cd_viscous = loads.cd_viscous;
			snapshot_.cs_viscous = loads.cs_viscous;
			snapshot_.cl_viscous = loads.cl_viscous;
			snapshot_.cd = loads.cd;
			snapshot_.cs = loads.cs;
			snapshot_.cl = loads.cl;
			snapshot_.conservation_valid = true;
			snapshot_.max_abs_divergence = conservation.max_abs_divergence;
			snapshot_.volume_weighted_rms_divergence = conservation.volume_weighted_rms_divergence;
			snapshot_.max_integrated_flux_error = conservation.max_integrated_flux_error;
			snapshot_.absolute_integrated_flux_error = conservation.absolute_integrated_flux_error;
			snapshot_.net_integrated_flux_error = conservation.net_integrated_flux_error;
			++snapshot_.surface_generation;
		}
		++snapshot_.generation;
	}

	void ParagliderSimWorker::run()
	{
		try
		{
			{
				std::lock_guard lock(recording_mutex_);recording_started_.store(recording_options_.any());
			}
			paracfd::core::ExternalAeroStepStats stats = core_->initialize();
			publishFlowField();
			// The zero-time projection creates a divergence-free, impermeable initial
			// velocity. Its pressure is the impulse required to create that state from
			// uniform flow, not an evolved aerodynamic pressure field.
			publish(stats, false);
			bool record_every_step=false;
			{
				std::lock_guard lock(recording_mutex_);record_every_step=recording_options_.any()&&
					recording_options_.every_step;
			}
			if(record_every_step)captureRecordingFrame(stats.physical_time,steps_);
			else captureScheduledRecordingFrames(stats.physical_time,steps_);
			if (!stats.pressure.converged)
			{
				char message[192];
				std::snprintf(message, sizeof(message),
					"initial pressure projection did not converge after %d iterations "
					"(relative residual %.6g; requested tolerance %.6g)",
					stats.pressure.iterations, stats.pressure.relative_residual,
					core_->config().solver.projection_tolerance);
				std::fprintf(stderr, "[paraglider-worker] ERROR: %s\n", message);
				std::fflush(stderr);
				throw std::runtime_error(message);
			}
			while (!stop_.load())
			{
				bool advance = playing_.load();
				const int requests = step_requests_.load();
				if (!advance && requests > 0)
				{
					advance = true;
					step_requests_.fetch_sub(1);
				}
				if (!advance)
				{
					QThread::msleep(10);
					continue;
				}
				stats = core_->step();
				++steps_;
				if(!stats.pressure.converged)
				{
					// This is a numerical failure, not a convergence endpoint or a manual
					// pause. Publish its solver diagnostics, then let the common error path
					// stop both single runs and AoA sweeps truthfully.
					publish(stats,false);
					char message[192];
					std::snprintf(message,sizeof(message),
						"pressure projection failed at evolved step %lld after %d iterations "
						"(relative residual %.6g; requested tolerance %.6g)",
						steps_,stats.pressure.iterations,stats.pressure.relative_residual,
						core_->config().solver.projection_tolerance);
					throw std::runtime_error(message);
				}
				// Surface loads require a deliberately throttled pressure download; scalar
				// solver telemetry remains available after every GPU step.
				bool record_due=false,record_every_step=false;
				{
					std::lock_guard lock(recording_mutex_);record_every_step=recording_options_.any()&&
						recording_options_.every_step;
					record_due=recording_options_.any()&&(record_every_step||
						next_recording_time_<=stats.physical_time+1e-10);
				}
				const bool convergence_publish=steps_==1||(steps_%10)==0;
				const bool publish_fields=convergence_publish||record_due;
				if (publish_fields) publishFlowField();
				publish(stats,publish_fields,convergence_publish);
				if(record_due)
				{
					if(record_every_step)captureRecordingFrame(stats.physical_time,steps_);
					else captureScheduledRecordingFrames(stats.physical_time,steps_);
				}
			}
		}
		catch (const std::exception& exception)
		{
			playing_.store(false);
			std::lock_guard lock(snapshot_mutex_);
			snapshot_.error = exception.what();
			snapshot_.converged = false;
			snapshot_.playing = false;
			++snapshot_.generation;
		}
		recording_started_.store(false);
		emit finished();
	}
}
