#pragma once

#include "core/geometry/tri_mesh.h"
#include "core/paraglider_config.h"

#include <QMainWindow>

#include <cstdint>
#include <memory>
#include <vector>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QMenu;
class QPushButton;
class QSlider;
class QSpinBox;
class QThread;
class QTimer;

namespace paracfd::core { class ExternalAeroCore; }

namespace paracfd::gui
{
	class ParagliderSimWorker;
	class SliceViewer;

	// Purpose-built paraglider application shell. It owns no channel, seabed, solid
	// voxel, or building state: STEP -> placed two-sided TriMesh -> static AMR/EB -> GPU.
	class ParagliderWindow final : public QMainWindow
	{
		Q_OBJECT
	public:
		explicit ParagliderWindow(QWidget* parent=nullptr);
		~ParagliderWindow() override;

		bool loadStepFile(const QString& path,bool infer_orientation=true,bool remember_file=true);
		bool loadConfigFile(const QString& path,bool build_after_load=false);
		bool saveConfigFile(const QString& path);
		bool loadLastFile();
		bool buildGrid();
		void setThinYDiagnostic(double span_fraction,double width_metres=0.125);
		SliceViewer* viewer() const{return viewer_;}
		long long steps()const{return last_steps_;}
		double physicalTime()const{return last_time_;}

	protected:
		void closeEvent(QCloseEvent* event)override;

	private:
		void buildMenus();
		void buildControls();
		void rememberRecentFile(const QString& path);
		void refreshRecentFiles();
		void configToUi(const paracfd::core::ParagliderConfig& config);
		paracfd::core::ParagliderConfig configFromUi()const;
		void normalizePlacementToDomain();
		void rotateWing(double degrees,const paracfd::core::Vec3d& axis);
		void updateGridReadout();
		void updateSnapshot();
		void applySurfaceColour();
		void updateDebugBoxes();
		void restoreFullWingDisplay();
		void spawnWorker(std::unique_ptr<paracfd::core::ExternalAeroCore> core);
		void shutdownWorker();

		SliceViewer* viewer_=nullptr;
		ParagliderSimWorker* worker_=nullptr;
		QThread* worker_thread_=nullptr;
		QTimer* repaint_timer_=nullptr;
		QMenu* recent_files_menu_=nullptr;

		paracfd::core::ParagliderConfig config_;
		paracfd::core::TriMesh source_mesh_;
		QString step_path_,config_path_;
		std::uint64_t snapshot_generation_=0,surface_generation_=0;
		long long last_steps_=0;
		double last_time_=0;

		QPushButton *build_button_=nullptr,*play_button_=nullptr,*step_button_=nullptr;
		QDoubleSpinBox *speed_=nullptr,*rho_=nullptr,*nu_=nullptr,*upstream_=nullptr,*downstream_=nullptr,
			*lateral_=nullptr,*vertical_=nullptr,*base_h_=nullptr,*wing_refine_=nullptr,*surface_refine_=nullptr,
			*wake_length_=nullptr,*wake_radius_=nullptr,*min_volume_fraction_=nullptr,*min_aperture_area_fraction_=nullptr,
			*cfl_=nullptr,*smagorinsky_=nullptr,*projection_tolerance_=nullptr,
			*reference_area_=nullptr,*reference_length_=nullptr,*tessellation_=nullptr;
		QSpinBox *levels_=nullptr,*brick_size_=nullptr,*projection_iterations_=nullptr;
		QComboBox *field_=nullptr,*slice_axis_=nullptr,*surface_colour_=nullptr,*arrow_mode_=nullptr,*tracer_mode_=nullptr;
		QSlider *slice_position_=nullptr,*auto_pause_sensitivity_=nullptr;
		QCheckBox *auto_range_=nullptr,*show_slice_=nullptr,*show_model_=nullptr,*show_amr_=nullptr,
			*show_eb_=nullptr,*show_arrows_=nullptr,*show_tracers_=nullptr,*thin_y_debug_=nullptr,*auto_pause_=nullptr;
		QDoubleSpinBox *thin_y_fraction_=nullptr,*thin_y_width_=nullptr;
		QLabel *wing_label_=nullptr,*grid_readout_=nullptr,*solver_readout_=nullptr,*load_readout_=nullptr;

		std::vector<std::array<float,6>> amr_boxes_,eb_boxes_;
		std::vector<float> cp_plus_,cp_minus_,delta_cp_;
		float delta_cp_range_=1,side_cp_range_=1;
		bool thin_debug_display_=false;
	};
}
