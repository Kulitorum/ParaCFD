#include "gui/paraglider_window.h"

#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/fluid/amr_pressure.h"
#include "core/fluid/external_aero_core.h"
#include "core/geometry/embedded_boundary.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "gui/paraglider_sim_worker.h"
#include "gui/display_info.h"
#include "gui/slice_viewer.h"

#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QPushButton>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSettings>
#include <QSlider>
#include <QSpinBox>
#include <QStatusBar>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>

namespace paracfd::gui
{
	using namespace paracfd::core;

	namespace
	{
		class ScrollSafeDoubleSpinBox final:public QDoubleSpinBox
		{
		protected:
			void wheelEvent(QWheelEvent* event)override{event->ignore();}
		};

		class ScrollSafeSpinBox final:public QSpinBox
		{
		protected:
			void wheelEvent(QWheelEvent* event)override{event->ignore();}
		};

		class ScrollSafeComboBox final:public QComboBox
		{
		protected:
			void wheelEvent(QWheelEvent* event)override{event->ignore();}
		};

		class ScrollSafeSlider final:public QSlider
		{
		public:
			using QSlider::QSlider;
		protected:
			void wheelEvent(QWheelEvent* event)override{event->ignore();}
		};

		QDoubleSpinBox* real_spin(double lo,double hi,double value,int decimals,const char* suffix="")
		{
			auto* spin=new ScrollSafeDoubleSpinBox;spin->setRange(lo,hi);spin->setDecimals(decimals);spin->setValue(value);
			spin->setSuffix(QString::fromUtf8(suffix));spin->setKeyboardTracking(false);return spin;
		}

		ModelPlacement left_rotation(const ModelPlacement& source,double degrees,const Vec3d& axis)
		{
			const double radians=degrees*3.14159265358979323846/180.0,c=std::cos(radians),s=std::sin(radians);
			const double length=std::sqrt(axis.x*axis.x+axis.y*axis.y+axis.z*axis.z);
			const double x=axis.x/length,y=axis.y/length,z=axis.z/length,t=1-c;
			const double r[9]={t*x*x+c,t*x*y-s*z,t*x*z+s*y,t*x*y+s*z,t*y*y+c,t*y*z-s*x,t*x*z-s*y,t*y*z+s*x,t*z*z+c};
			ModelPlacement out=source;
			for(int row=0;row<3;++row)for(int column=0;column<3;++column)
				out.m[3*row+column]=r[3*row+0]*source.m[column]+r[3*row+1]*source.m[3+column]+r[3*row+2]*source.m[6+column];
			out.tx=out.ty=out.tz=0;return out;
		}

		QString resolve_step_path(const QString& config_path,const std::string& stored)
		{
			const QString candidate=QString::fromStdString(stored);
			const QFileInfo candidate_info(candidate);
			if(candidate_info.isAbsolute())return candidate_info.absoluteFilePath();

			// Relative STEP provenance is primarily relative to the config.  Historic
			// ParaCFD configs used repository-root paths such as Test-Data/wing.step,
			// though, so also search the config directory's ancestors.  This keeps a
			// pinned Release executable independent of its process working directory.
			QDir directory=QFileInfo(config_path).absoluteDir();
			for(;;)
			{
				const QString resolved=directory.absoluteFilePath(candidate);
				if(QFileInfo::exists(resolved))return QFileInfo(resolved).absoluteFilePath();
				if(!directory.cdUp())break;
			}
			return QFileInfo(config_path).absoluteDir().absoluteFilePath(candidate);
		}
	}

	ParagliderWindow::ParagliderWindow(QWidget* parent):QMainWindow(parent)
	{
		setWindowTitle("ParaCFD — GPU paraglider aerodynamics");resize(1320,820);
		viewer_=new SliceViewer(this);setCentralWidget(viewer_);
		SimInfo initial;initial.nx=initial.ny=initial.nz=32;initial.h=0.5;initial.coarse_h=0.5;initial.finest_h=0.5;initial.Lx=initial.Ly=initial.Lz=16;initial.U=config_.freestream.speed;initial.rho=config_.freestream.rho;initial.nu=config_.freestream.nu;initial.name="paraglider";viewer_->setInfo(initial);
		viewer_->setShowSlice(true);viewer_->setShowModel(true);viewer_->setShowArrows(true);viewer_->setArrowMode3D(true);viewer_->setArrowDensity(400);viewer_->setArrowSpeedMult(0.1f);viewer_->setArrowSizeMult(0.5f);viewer_->setShowTracers(true);viewer_->setTracerMode3D(true);viewer_->setTracerGridDensity(7);
		buildMenus();buildControls();configToUi(config_);
		repaint_timer_=new QTimer(this);connect(repaint_timer_,&QTimer::timeout,this,[this]{updateSnapshot();viewer_->update();});repaint_timer_->start(16);
		statusBar()->showMessage("Open a STEP wing, inspect its orientation, then Build CFD Grid.");
	}

	ParagliderWindow::~ParagliderWindow(){shutdownWorker();}

	void ParagliderWindow::buildMenus()
	{
		auto* file=menuBar()->addMenu("&File");
		auto* open_step=file->addAction("Open STEP...");connect(open_step,&QAction::triggered,this,[this]{const QString path=QFileDialog::getOpenFileName(this,"Open paraglider STEP",{},"STEP files (*.step *.stp);;All files (*)");if(!path.isEmpty())loadStepFile(path,true);});
		auto* open_config=file->addAction("Open paraglider config...");connect(open_config,&QAction::triggered,this,[this]{const QString path=QFileDialog::getOpenFileName(this,"Open ParaCFD config",{},"JSON files (*.json);;All files (*)");if(!path.isEmpty())loadConfigFile(path,false);});
		auto* save_config=file->addAction("Save paraglider config as...");connect(save_config,&QAction::triggered,this,[this]{QString path=QFileDialog::getSaveFileName(this,"Save ParaCFD config",config_path_,"JSON files (*.json)");if(!path.isEmpty()){if(!path.endsWith(".json",Qt::CaseInsensitive))path+=".json";saveConfigFile(path);}});
		recent_files_menu_=file->addMenu("Recent files");refreshRecentFiles();
		file->addSeparator();auto* quit=file->addAction("Exit");connect(quit,&QAction::triggered,this,&QWidget::close);
	}

	void ParagliderWindow::rememberRecentFile(const QString& path)
	{
		const QFileInfo info(path);const QString absolute=info.canonicalFilePath().isEmpty()?info.absoluteFilePath():info.canonicalFilePath();if(absolute.isEmpty())return;QSettings settings;QStringList recent=settings.value("recentFiles").toStringList();recent.removeIf([&](const QString& entry){return QString::compare(QFileInfo(entry).absoluteFilePath(),absolute,Qt::CaseInsensitive)==0;});recent.prepend(absolute);while(recent.size()>10)recent.removeLast();settings.setValue("recentFiles",recent);refreshRecentFiles();
	}

	void ParagliderWindow::refreshRecentFiles()
	{
		if(!recent_files_menu_)return;recent_files_menu_->clear();QSettings settings;QStringList recent=settings.value("recentFiles").toStringList(),valid;for(const QString& path:recent)if(QFileInfo::exists(path)&&!valid.contains(path,Qt::CaseInsensitive))valid.push_back(QFileInfo(path).absoluteFilePath());if(valid!=recent)settings.setValue("recentFiles",valid);if(valid.empty()){auto* empty=recent_files_menu_->addAction("No recent files");empty->setEnabled(false);}else for(int index=0;index<valid.size();++index){const QFileInfo info(valid[index]);auto* action=recent_files_menu_->addAction(QString("&%1  %2 — %3").arg(index+1).arg(info.fileName(),QDir::toNativeSeparators(info.absolutePath())));action->setToolTip(QDir::toNativeSeparators(valid[index]));connect(action,&QAction::triggered,this,[this,path=valid[index]]{if(path.endsWith(".json",Qt::CaseInsensitive))loadConfigFile(path,false);else loadStepFile(path,true);});}recent_files_menu_->addSeparator();auto* clear=recent_files_menu_->addAction("Clear recent files");clear->setEnabled(!valid.empty());connect(clear,&QAction::triggered,this,[this]{QSettings{}.remove("recentFiles");refreshRecentFiles();});
	}

	bool ParagliderWindow::loadLastFile()
	{
		QSettings settings;const QStringList recent=settings.value("recentFiles").toStringList();for(const QString& path:recent)if(QFileInfo::exists(path)){const bool loaded=path.endsWith(".json",Qt::CaseInsensitive)?loadConfigFile(path,false):loadStepFile(path,true);if(loaded){statusBar()->showMessage(QString("Restored recent wing: %1").arg(QFileInfo(path).fileName()),6000);return true;}}refreshRecentFiles();return false;
	}

	void ParagliderWindow::buildControls()
	{
		auto* dock=new QDockWidget("Paraglider CFD",this);dock->setAllowedAreas(Qt::LeftDockWidgetArea|Qt::RightDockWidgetArea);dock->setFeatures(QDockWidget::DockWidgetMovable|QDockWidget::DockWidgetFloatable);
		auto* panel=new QWidget;auto* column=new QVBoxLayout(panel);column->setContentsMargins(8,8,8,8);
		wing_label_=new QLabel("No STEP wing loaded");wing_label_->setWordWrap(true);column->addWidget(wing_label_);
		auto* wind_hint=new QLabel(QString::fromUtf8("FREESTREAM +X  →\nBBox infers span/chord; leading/trailing needs confirmation."));wind_hint->setStyleSheet("font-weight:600;color:#ff9a75;background:#20242a;padding:6px;");column->addWidget(wind_hint);
		auto* orientation=new QGridLayout;auto* flip=new QPushButton("Flip LE/TE 180°");auto* yaw=new QPushButton("Yaw +90°");auto* aoa_up=new QPushButton("AoA +1°");auto* aoa_down=new QPushButton("AoA -1°");connect(flip,&QPushButton::clicked,this,[this]{rotateWing(180,{0,0,1});});connect(yaw,&QPushButton::clicked,this,[this]{rotateWing(90,{0,0,1});});connect(aoa_up,&QPushButton::clicked,this,[this]{rotateWing(1,{0,1,0});});connect(aoa_down,&QPushButton::clicked,this,[this]{rotateWing(-1,{0,1,0});});orientation->addWidget(flip,0,0);orientation->addWidget(yaw,0,1);orientation->addWidget(aoa_up,1,0);orientation->addWidget(aoa_down,1,1);column->addLayout(orientation);

		auto* run_group=new QGroupBox("Run");auto* run_column=new QVBoxLayout(run_group);build_button_=new QPushButton("Build CFD Grid + Run");build_button_->setToolTip("Build the static AMR/embedded-boundary grid, initialize the pressure field, and immediately run the CFD simulation.");play_button_=new QPushButton("Play");play_button_->setCheckable(true);step_button_=new QPushButton("Step");play_button_->setEnabled(false);step_button_->setEnabled(false);connect(build_button_,&QPushButton::clicked,this,[this]{buildGrid();});connect(play_button_,&QPushButton::toggled,this,[this](bool on){play_button_->setText(on?"Pause":"Play");if(worker_)worker_->setPlaying(on);});connect(step_button_,&QPushButton::clicked,this,[this]{if(worker_)worker_->stepOnce();});auto* run_row=new QHBoxLayout;run_row->addWidget(play_button_);run_row->addWidget(step_button_);run_column->addWidget(build_button_);run_column->addLayout(run_row);column->addWidget(run_group);
		auto* thin_group=new QGroupBox("Cropped Y-span diagnosis");auto* thin_form=new QFormLayout(thin_group);thin_y_debug_=new QCheckBox("Run cropped Y volume");thin_y_debug_->setToolTip("Crop the oriented wing and CFD domain to a configurable Y width at one span station. This uses the real uniform EB solver but is a diagnostic case, not a physical whole-wing result.");thin_y_fraction_=real_spin(0.05,0.95,0.5,3);thin_y_fraction_->setSingleStep(0.025);thin_y_fraction_->setToolTip("Span station through the oriented wing bbox: 0 is one tip, 0.5 is centre, 1 is the other tip.");thin_y_width_=real_spin(0.001,100.0,0.125,4," m");thin_y_width_->setSingleStep(0.125);thin_y_width_->setToolTip("Requested physical Y width. It is rounded upward to an integer number of configured finest cells; the actual width is shown below.");thin_form->addRow(thin_y_debug_);thin_form->addRow("Span station",thin_y_fraction_);thin_form->addRow("Y volume width",thin_y_width_);column->addWidget(thin_group);connect(thin_y_debug_,&QCheckBox::toggled,this,[this]{updateGridReadout();});connect(thin_y_fraction_,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,[this]{updateGridReadout();});connect(thin_y_width_,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,[this]{updateGridReadout();});

		auto* physics=new QGroupBox("Freestream / solver (Build resets)");auto* form=new QFormLayout(physics);
		speed_=real_spin(0,100,10,3," m/s");rho_=real_spin(0.1,10,1.225,4," kg/m³");nu_=real_spin(1e-8,1e-2,1.5e-5,8," m²/s");tessellation_=real_spin(0.01,100,2,2," mm");
		upstream_=real_spin(0,100,3,2," m");downstream_=real_spin(0,200,8,2," m");lateral_=real_spin(0,100,3,2," m");vertical_=real_spin(0,100,3,2," m");
		base_h_=real_spin(0.015625,10,0.25,5," m");levels_=new ScrollSafeSpinBox;levels_->setRange(1,6);brick_size_=new ScrollSafeSpinBox;brick_size_->setRange(8,64);brick_size_->setSingleStep(8);
		wing_refine_=real_spin(0,50,1,3," m");surface_refine_=real_spin(0,20,0.35,3," m");wake_length_=real_spin(0,200,8,2," m");wake_radius_=real_spin(0,100,2,2," m");
		min_volume_fraction_=real_spin(0.01,0.49,0.25,3);min_volume_fraction_->setSingleStep(0.01);min_aperture_area_fraction_=real_spin(0,0.1,1e-4,6);min_aperture_area_fraction_->setSingleStep(1e-4);
		cfl_=real_spin(0.02,0.95,0.7,2);smagorinsky_=real_spin(0,0.4,0.1,3);projection_tolerance_=real_spin(1e-8,1e-2,1e-5,8);projection_iterations_=new ScrollSafeSpinBox;projection_iterations_->setRange(20,5000);projection_iterations_->setSingleStep(50);
		reference_area_=real_spin(0,10000,0,3," m²");reference_length_=real_spin(0,1000,0,3," m");
		form->addRow("Speed",speed_);form->addRow("Density",rho_);form->addRow("Kinematic viscosity",nu_);form->addRow("STEP deflection",tessellation_);form->addRow("Upstream margin",upstream_);form->addRow("Downstream margin",downstream_);form->addRow("Lateral margin",lateral_);form->addRow("Vertical margin",vertical_);form->addRow("Base cell size",base_h_);form->addRow("AMR levels",levels_);form->addRow("Brick size",brick_size_);form->addRow("Wing refine distance",wing_refine_);form->addRow("Surface refine distance",surface_refine_);form->addRow("Wake length",wake_length_);form->addRow("Wake radius",wake_radius_);form->addRow("Min fragment volume / h³",min_volume_fraction_);form->addRow("Min aperture area / h²",min_aperture_area_fraction_);form->addRow("CFL",cfl_);form->addRow("Smagorinsky Cs",smagorinsky_);form->addRow("Projection tolerance",projection_tolerance_);form->addRow("Projection max iterations",projection_iterations_);form->addRow("Reference area",reference_area_);form->addRow("Reference length",reference_length_);column->addWidget(physics);
		for(auto* spin:{speed_,rho_,nu_,tessellation_,upstream_,downstream_,lateral_,vertical_,base_h_,wing_refine_,surface_refine_,wake_length_,wake_radius_,min_volume_fraction_,min_aperture_area_fraction_,cfl_,smagorinsky_,projection_tolerance_,reference_area_,reference_length_})connect(spin,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,[this]{updateGridReadout();});connect(levels_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this]{updateGridReadout();});connect(brick_size_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this]{updateGridReadout();});
		grid_readout_=new QLabel;grid_readout_->setStyleSheet("font-family:Consolas;color:#000;");grid_readout_->setWordWrap(true);column->addWidget(grid_readout_);

		auto* visualization=new QGroupBox("Visualization");auto* viz_form=new QFormLayout(visualization);field_=new ScrollSafeComboBox;field_->addItems({"Speed |u|","X velocity u","Y velocity v","Z velocity w","Pressure p"});slice_axis_=new ScrollSafeComboBox;slice_axis_->addItems({"X-normal","Y-normal","Z-normal"});slice_axis_->setCurrentIndex(2);slice_position_=new ScrollSafeSlider(Qt::Horizontal);slice_position_->setRange(0,1000);slice_position_->setValue(500);auto_range_=new QCheckBox("Auto range");auto_range_->setChecked(true);show_slice_=new QCheckBox("Slice");show_slice_->setChecked(true);show_model_=new QCheckBox("STEP surface");show_model_->setChecked(true);show_amr_=new QCheckBox("AMR grid");show_amr_->setToolTip("Show or hide the cyan three-dimensional AMR brick boundaries.");show_amr_->setChecked(false);show_eb_=new QCheckBox("EB cells");show_eb_->setChecked(false);show_arrows_=new QCheckBox("Velocity arrows");show_arrows_->setChecked(true);show_tracers_=new QCheckBox("Flow tracers");show_tracers_->setChecked(true);surface_colour_=new ScrollSafeComboBox;surface_colour_->addItems({"Visible side Cp+/Cp-",QString::fromUtf8("Pressure difference ΔCp"),"Plus side Cp+","Minus side Cp-"});surface_colour_->setCurrentIndex(1);auto* layers=new QWidget;auto* layers_grid=new QGridLayout(layers);layers_grid->setContentsMargins(0,0,0,0);const QList<QCheckBox*> layer_checks={show_slice_,show_model_,show_amr_,show_eb_,show_arrows_,show_tracers_};for(int q=0;q<layer_checks.size();++q)layers_grid->addWidget(layer_checks[q],q/2,q%2);auto* views=new QWidget;auto* views_column=new QVBoxLayout(views);views_column->setContentsMargins(0,0,0,0);auto* fit_wing=new QPushButton("Fit Wing");auto* fit_domain=new QPushButton("Fit Domain");views_column->addWidget(fit_wing);views_column->addWidget(fit_domain);viz_form->addRow("Field",field_);viz_form->addRow("Slice plane",slice_axis_);viz_form->addRow("Plane position",slice_position_);viz_form->addRow(auto_range_);viz_form->addRow("View",views);viz_form->addRow("Layers",layers);viz_form->addRow("Canopy colour",surface_colour_);column->addWidget(visualization);
		connect(field_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int index){static constexpr Field fields[]={Field::SpeedMag,Field::VelU,Field::VelV,Field::VelW,Field::Pressure};viewer_->setField(fields[index]);});connect(slice_axis_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int index){viewer_->setAxis(static_cast<Axis>(index));});connect(slice_position_,&QSlider::valueChanged,this,[this](int value){viewer_->setPlaneFraction(value/1000.0f);});connect(auto_range_,&QCheckBox::toggled,viewer_,&SliceViewer::setAutoRange);connect(fit_wing,&QPushButton::clicked,viewer_,&SliceViewer::frameWingView);connect(fit_domain,&QPushButton::clicked,viewer_,&SliceViewer::frameDomainView);connect(show_slice_,&QCheckBox::toggled,viewer_,&SliceViewer::setShowSlice);connect(show_model_,&QCheckBox::toggled,viewer_,&SliceViewer::setShowModel);connect(show_amr_,&QCheckBox::toggled,this,[this]{updateDebugBoxes();});connect(show_eb_,&QCheckBox::toggled,this,[this]{updateDebugBoxes();});connect(show_arrows_,&QCheckBox::toggled,viewer_,&SliceViewer::setShowArrows);connect(show_tracers_,&QCheckBox::toggled,viewer_,&SliceViewer::setShowTracers);connect(surface_colour_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this]{applySurfaceColour();});

		auto* arrows_group=new QGroupBox("Flow arrows (visual only)");auto* arrows_form=new QFormLayout(arrows_group);arrow_mode_=new ScrollSafeComboBox;arrow_mode_->addItems({"3D volume","2D slice plane"});arrow_mode_->setCurrentIndex(viewer_->arrowMode3D()?0:1);arrow_mode_->setToolTip("3D advects arrows through the volume; 2D confines them to the selected slice and draws them as a depth-independent overlay.");auto* arrow_density=new ScrollSafeSlider(Qt::Horizontal);arrow_density->setRange(100,8000);arrow_density->setValue(viewer_->arrowDensity());arrow_density->setToolTip("Number of animated flow arrows.");auto* arrow_speed=new ScrollSafeSlider(Qt::Horizontal);arrow_speed->setRange(0,100);arrow_speed->setValue(static_cast<int>(std::lround(viewer_->arrowSpeedMult()*100.0f)));arrow_speed->setToolTip("Animation speed only: 0 freezes the arrows and 100 is full visual speed. This does not alter the CFD solution.");auto* arrow_size=new ScrollSafeSlider(Qt::Horizontal);arrow_size->setRange(5,100);arrow_size->setValue(static_cast<int>(std::lround(viewer_->arrowSizeMult()*50.0f)));arrow_size->setToolTip("Uniform arrow size: scales both length and width.");arrows_form->addRow("Mode",arrow_mode_);arrows_form->addRow("Density",arrow_density);arrows_form->addRow("Animation speed",arrow_speed);arrows_form->addRow("Size",arrow_size);column->addWidget(arrows_group);
		connect(arrow_mode_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int value){viewer_->setArrowMode3D(value==0);});connect(arrow_density,&QSlider::valueChanged,viewer_,&SliceViewer::setArrowDensity);connect(arrow_speed,&QSlider::valueChanged,this,[this](int value){viewer_->setArrowSpeedMult(value/100.0f);});connect(arrow_size,&QSlider::valueChanged,this,[this](int value){viewer_->setArrowSizeMult(value/50.0f);});

		auto* tracers_group=new QGroupBox("Flow tracers (visual only)");auto* tracers_form=new QFormLayout(tracers_group);tracer_mode_=new ScrollSafeComboBox;tracer_mode_->addItems({"3D volume","2D slice plane"});tracer_mode_->setCurrentIndex(viewer_->tracerMode3D()?0:1);tracer_mode_->setToolTip("3D seeds the inlet plane and integrates full 3D streamlines; 2D stays on the selected slice.");auto* tracer_density=new ScrollSafeSlider(Qt::Horizontal);tracer_density->setRange(2,120);tracer_density->setValue(viewer_->tracerGridDensity());tracer_density->setToolTip("Inlet seed count along the larger inlet dimension.");auto* tracer_length=new ScrollSafeSlider(Qt::Horizontal);tracer_length->setRange(50,2000);tracer_length->setValue(viewer_->tracerTrail());tracer_length->setToolTip("Maximum streamline length in integration steps.");auto* tracer_width=new ScrollSafeSlider(Qt::Horizontal);tracer_width->setRange(10,60);tracer_width->setValue(static_cast<int>(std::lround(viewer_->tracerWidth()*10.0f)));tracer_width->setToolTip("Tracer ribbon thickness in screen pixels.");auto* tracer_boring=new ScrollSafeSlider(Qt::Horizontal);tracer_boring->setRange(0,100000);tracer_boring->setSingleStep(1);tracer_boring->setPageStep(1000);tracer_boring->setTracking(true);tracer_boring->setValue(static_cast<int>(std::lround(100000.0f*(viewer_->tracerBoring()-0.99f)/0.03f)));tracer_boring->setToolTip("Hide paths by straightness (arc length / endpoint distance). Continuous linear range: 0.990 to 1.020; straight is 1.000.");tracers_form->addRow("Mode",tracer_mode_);tracers_form->addRow("Seed density",tracer_density);tracers_form->addRow("Length",tracer_length);tracers_form->addRow("Width",tracer_width);tracers_form->addRow("Hide boring",tracer_boring);column->addWidget(tracers_group);
		connect(tracer_mode_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int value){viewer_->setTracerMode3D(value==0);});connect(tracer_density,&QSlider::valueChanged,viewer_,&SliceViewer::setTracerGridDensity);connect(tracer_length,&QSlider::valueChanged,viewer_,&SliceViewer::setTracerTrail);connect(tracer_width,&QSlider::valueChanged,this,[this](int value){viewer_->setTracerWidth(value/10.0f);});connect(tracer_boring,&QSlider::valueChanged,this,[this](int value){viewer_->setTracerBoring(0.99f+0.03f*(value/100000.0f));});connect(tracer_boring,&QSlider::sliderPressed,this,[this]{viewer_->setTracerBoringInstant(true);});connect(tracer_boring,&QSlider::sliderReleased,this,[this]{viewer_->setTracerBoringInstant(false);});

		solver_readout_=new QLabel("Grid not built");solver_readout_->setWordWrap(true);solver_readout_->setStyleSheet("font-family:Consolas;color:#000;");load_readout_=new QLabel("Pressure loads unavailable");load_readout_->setWordWrap(true);load_readout_->setStyleSheet("font-family:Consolas;color:#000;");column->addWidget(solver_readout_);column->addWidget(load_readout_);column->addStretch();
		auto* scroll=new QScrollArea;scroll->setWidgetResizable(true);scroll->setWidget(panel);dock->setWidget(scroll);addDockWidget(Qt::LeftDockWidgetArea,dock);resizeDocks({dock},{380},Qt::Horizontal);
	}

	void ParagliderWindow::configToUi(const ParagliderConfig& c)
	{
		config_=c;speed_->setValue(c.freestream.speed);rho_->setValue(c.freestream.rho);nu_->setValue(c.freestream.nu);tessellation_->setValue(c.tessellation_deflection_mm);upstream_->setValue(c.domain.upstream_margin);downstream_->setValue(c.domain.downstream_margin);lateral_->setValue(c.domain.lateral_margin);vertical_->setValue(c.domain.vertical_margin);base_h_->setValue(c.amr.base_cell_size);levels_->setValue(c.amr.max_levels);brick_size_->setValue(c.amr.brick_size);wing_refine_->setValue(c.amr.wing_refinement_distance);surface_refine_->setValue(c.amr.surface_refinement_distance);wake_length_->setValue(c.amr.wake_length);wake_radius_->setValue(c.amr.wake_radius);min_volume_fraction_->setValue(c.amr.min_volume_fraction);min_aperture_area_fraction_->setValue(c.amr.min_aperture_area_fraction);cfl_->setValue(c.solver.cfl);smagorinsky_->setValue(c.solver.smagorinsky_cs);projection_tolerance_->setValue(c.solver.projection_tolerance);projection_iterations_->setValue(c.solver.projection_max_iterations);reference_area_->setValue(c.reference.area);reference_length_->setValue(c.reference.length);viewer_->setReferenceU(c.freestream.speed);updateGridReadout();
	}

	ParagliderConfig ParagliderWindow::configFromUi()const
	{
		ParagliderConfig c=config_;c.step_path=step_path_.toStdString();c.placement=config_.placement;c.tessellation_deflection_mm=tessellation_->value();c.freestream={speed_->value(),rho_->value(),nu_->value()};c.domain={upstream_->value(),downstream_->value(),lateral_->value(),vertical_->value()};c.amr.base_cell_size=base_h_->value();c.amr.max_levels=levels_->value();c.amr.brick_size=brick_size_->value();c.amr.wing_refinement_distance=wing_refine_->value();c.amr.surface_refinement_distance=surface_refine_->value();c.amr.wake_length=wake_length_->value();c.amr.wake_radius=wake_radius_->value();c.amr.min_volume_fraction=min_volume_fraction_->value();c.amr.min_aperture_area_fraction=min_aperture_area_fraction_->value();c.solver.cfl=cfl_->value();c.solver.smagorinsky_cs=smagorinsky_->value();c.solver.projection_tolerance=projection_tolerance_->value();c.solver.projection_max_iterations=projection_iterations_->value();c.reference.area=reference_area_->value();c.reference.length=reference_length_->value();return c;
	}

	void ParagliderWindow::setThinYDiagnostic(double span_fraction,double width_metres)
	{
		if(!thin_y_debug_||!thin_y_fraction_||!thin_y_width_)return;thin_y_fraction_->setValue(std::clamp(span_fraction,0.05,0.95));thin_y_width_->setValue(std::max(0.001,width_metres));thin_y_debug_->setChecked(true);
	}

	void ParagliderWindow::normalizePlacementToDomain()
	{
		if(source_mesh_.empty())return;const ParagliderConfig c=configFromUi();config_.placement=frame_wing_for_external_domain(source_mesh_,config_.placement,c.domain.upstream_margin,c.domain.lateral_margin,c.domain.vertical_margin,c.amr.base_cell_size*c.amr.brick_size);viewer_->setMeshPlacement(config_.placement);
	}

	void ParagliderWindow::rotateWing(double degrees,const Vec3d& axis)
	{
		if(source_mesh_.empty())return;shutdownWorker();restoreFullWingDisplay();config_.placement=left_rotation(config_.placement,degrees,axis);normalizePlacementToDomain();play_button_->setChecked(false);play_button_->setEnabled(false);step_button_->setEnabled(false);amr_boxes_.clear();eb_boxes_.clear();updateDebugBoxes();viewer_->clearTriangleSurfaceColouring();statusBar()->showMessage(QString("Wing rotated %1°; rebuild the static CFD grid.").arg(degrees),5000);
	}

	bool ParagliderWindow::loadStepFile(const QString& path,bool infer_orientation,bool remember_file)
	{
		shutdownWorker();
		QApplication::setOverrideCursor(Qt::WaitCursor);
		std::string error;
		TriMesh mesh=load_step_mesh(path.toStdString(),tessellation_->value(),&error);
		QApplication::restoreOverrideCursor();
		if(mesh.empty())
		{
			QMessageBox::critical(this,"STEP import failed",QString::fromStdString(error));
			return false;
		}

		source_mesh_=mesh;
		thin_debug_display_=false;viewer_->setThinDebugState(false,0);viewer_->setSimulationStep(0);
		step_path_=QFileInfo(path).absoluteFilePath();
		config_.step_path=step_path_.toStdString();
		QString orientation_note="placement restored from config";
		if(infer_orientation)
		{
			config_.placement=ModelPlacement{};
			const double dx=mesh.bbox_max[0]-mesh.bbox_min[0];
			const double dy=mesh.bbox_max[1]-mesh.bbox_min[1];
			const HorizontalWingAxes axes=infer_horizontal_wing_axes(mesh);
			if(axes.valid)
			{
				config_.placement=left_rotation(config_.placement,axes.yaw_degrees,{0,0,1});
				orientation_note=axes.span_axis==0
					?"bbox: span X / chord Y; assumed leading edge -Y -> upstream -X"
					:"bbox: span Y / chord X; assumed leading edge -X -> upstream -X";
			}
			else orientation_note="bbox span/chord ambiguous; imported axes retained";
			std::fprintf(stderr,"[paraglider-orientation] bbox %.3f x %.3f m: %s\n",
				dx,dy,orientation_note.toUtf8().constData());
		}

		viewer_->setMesh(source_mesh_);
		normalizePlacementToDomain();
		viewer_->frameWingView();
		wing_label_->setText(QString("%1\n%2 face triangles — wires ignored\n%3\nLE/TE polarity must be confirmed")
			.arg(QFileInfo(path).fileName()).arg(source_mesh_.triangle_count()).arg(orientation_note));
		amr_boxes_.clear();eb_boxes_.clear();cp_plus_.clear();cp_minus_.clear();delta_cp_.clear();
		updateDebugBoxes();viewer_->clearTriangleSurfaceColouring();
		play_button_->setChecked(false);play_button_->setEnabled(false);step_button_->setEnabled(false);
		if(remember_file)rememberRecentFile(step_path_);
		statusBar()->showMessage("STEP loaded and face bbox centred. Confirm leading/trailing direction, then Build CFD Grid.",8000);
		return true;
	}

	bool ParagliderWindow::loadConfigFile(const QString& path,bool build_after_load)
	{
		ParagliderConfig loaded;std::string error;if(!load_paraglider_config(path.toStdString(),loaded,&error)){QMessageBox::critical(this,"Config load failed",QString::fromStdString(error));return false;}config_path_=QFileInfo(path).absoluteFilePath();configToUi(loaded);if(loaded.step_path.empty()){statusBar()->showMessage("Config loaded; choose its STEP wing.",6000);return true;}const QString step=resolve_step_path(config_path_,loaded.step_path);if(!loadStepFile(step,false,false))return false;config_.placement=loaded.placement;normalizePlacementToDomain();viewer_->setMeshPlacement(config_.placement);viewer_->frameWingView();rememberRecentFile(config_path_);return !build_after_load||buildGrid();
	}

	bool ParagliderWindow::saveConfigFile(const QString& path)
	{
		config_=configFromUi();if(!thin_debug_display_)config_.placement=viewer_->modelPlacement();std::string error;if(!save_paraglider_config(path.toStdString(),config_,&error)){QMessageBox::critical(this,"Config save failed",QString::fromStdString(error));return false;}config_path_=QFileInfo(path).absoluteFilePath();rememberRecentFile(config_path_);statusBar()->showMessage(QString("Saved %1").arg(config_path_),5000);return true;
	}

	void ParagliderWindow::updateGridReadout()
	{
		if(!grid_readout_)return;const auto c=configFromUi();const int ratio=1<<std::max(0,c.amr.max_levels-1);const double finest=c.amr.base_cell_size/ratio;QString text=QString("finest h = %1 m (2:1 ×%2)\nbrick = %3³, coarse width %4 m\nEB merge < %5 h³, aperture cutoff < %6 h²\nproduction fields = %7").arg(finest,0,'g',5).arg(ratio).arg(c.amr.brick_size).arg(c.amr.base_cell_size*c.amr.brick_size,0,'g',5).arg(c.amr.min_volume_fraction,0,'g',4).arg(c.amr.min_aperture_area_fraction,0,'g',4).arg(sizeof(Real)==4?"FP32":"FP64 validation");if(thin_y_debug_&&thin_y_debug_->isChecked()){const int layers=std::max(2,static_cast<int>(std::ceil(thin_y_width_->value()/finest-1e-9)));text+=QString("\nCROPPED-Y DEBUG: uniform h=%1 m, %2 cells = %3 m, station=%4%5").arg(finest,0,'g',5).arg(layers).arg(layers*finest,0,'g',5).arg(thin_y_fraction_->value(),0,'f',3).arg(layers>128?" (too wide: max 128 cells)":"");}grid_readout_->setText(text);
	}

	void ParagliderWindow::restoreFullWingDisplay()
	{
		if(!thin_debug_display_||source_mesh_.empty())return;viewer_->setMesh(source_mesh_);viewer_->setMeshPlacement(config_.placement);viewer_->setThinDebugState(false,0);thin_debug_display_=false;
	}

	bool ParagliderWindow::buildGrid()
	{
		if(source_mesh_.empty()){QMessageBox::information(this,"No wing","Open a STEP wing first.");return false;}
		restoreFullWingDisplay();shutdownWorker();play_button_->setChecked(false);play_button_->setEnabled(false);step_button_->setEnabled(false);config_=configFromUi();config_.placement=viewer_->modelPlacement();normalizePlacementToDomain();config_.placement=viewer_->modelPlacement();
		const bool thin_debug=thin_y_debug_&&thin_y_debug_->isChecked();int debug_layers=0;double debug_width=0;ParagliderConfig run_config=config_;TriMesh wing;
		if(thin_debug)
		{
			const int ratio=1<<std::max(0,config_.amr.max_levels-1);const double h=config_.amr.base_cell_size/ratio;debug_layers=std::max(2,static_cast<int>(std::ceil(thin_y_width_->value()/h-1e-9)));if(debug_layers>128){QMessageBox::critical(this,"Cropped Y volume too wide",QString("The requested width needs %1 finest cells. This diagnostic supports at most 128; reduce the width or use the full-wing run.").arg(debug_layers));return false;}debug_width=debug_layers*h;ModelPlacement orientation=config_.placement;orientation.tx=orientation.ty=orientation.tz=0;const TriMesh oriented=placed_mesh(source_mesh_,orientation);const double fraction=thin_y_fraction_->value(),centre=oriented.bbox_min[1]+fraction*(oriented.bbox_max[1]-oriented.bbox_min[1]);TriMesh clipped=clip_mesh_to_axis_slab(oriented,1,centre-0.5*debug_width,centre+0.5*debug_width);if(clipped.empty()){QMessageBox::critical(this,"Cropped Y volume failed","The selected Y slab contains no fabric triangles.");return false;}run_config.amr.base_cell_size=h;run_config.amr.max_levels=1;run_config.amr.brick_size=debug_layers;run_config.domain.lateral_margin=0;run_config.reference.area=0;run_config.reference.length=0;const ModelPlacement frame=frame_wing_for_external_domain(clipped,ModelPlacement{},run_config.domain.upstream_margin,0,run_config.domain.vertical_margin,debug_width);wing=placed_mesh(clipped,frame);run_config.placement=ModelPlacement{};std::fprintf(stderr,"[paraglider-thin-y] station=%.6g source-y=%.6g width=%.6g m h=%.6g layers=%d triangles=%zu\n",fraction,centre,debug_width,h,debug_layers,wing.triangle_count());
		}
		else wing=placed_mesh(source_mesh_,config_.placement);
		std::fprintf(stderr,"[paraglider-placement] t=[%.17g %.17g %.17g] M=[%.17g %.17g %.17g; %.17g %.17g %.17g; %.17g %.17g %.17g]\n",config_.placement.tx,config_.placement.ty,config_.placement.tz,config_.placement.m[0],config_.placement.m[1],config_.placement.m[2],config_.placement.m[3],config_.placement.m[4],config_.placement.m[5],config_.placement.m[6],config_.placement.m[7],config_.placement.m[8]);TriangleBvh bvh(wing);
		QApplication::setOverrideCursor(Qt::WaitCursor);std::unique_ptr<ExternalAeroCore> core;try{core=std::make_unique<ExternalAeroCore>(wing,bvh,run_config);}catch(const std::exception& e){QApplication::restoreOverrideCursor();QMessageBox::critical(this,"CFD grid failed",e.what());return false;}QApplication::restoreOverrideCursor();const AmrHierarchy& hierarchy=core->hierarchy();
		amr_boxes_.clear();eb_boxes_.clear();for(const auto& level:hierarchy.levels())for(const auto& brick:level.bricks)if(brick.active()){const float width=hierarchy.brick_size()*brick.h;amr_boxes_.push_back({(float)brick.origin.x,(float)brick.origin.y,(float)brick.origin.z,(float)brick.origin.x+width,(float)brick.origin.y+width,(float)brick.origin.z+width});}
		std::size_t fragments=0,apertures=0,patches=0,unresolved=0,static_pockets=0,discarded_apertures=0;double discarded_aperture_area=0;for(const auto& level:core->embedded_boundary().levels){const auto& eb=level.topology;discarded_apertures+=eb.discarded_subgrid_apertures;discarded_aperture_area+=eb.discarded_subgrid_aperture_area;for(const auto& fragment:eb.fragments)if(level.owned_cell[fragment.parent_cell]){++fragments;if(fragment.pressure_static)++static_pockets;}for(const auto& aperture:eb.apertures)if(level.owned_cell[aperture.parent_face_cell])++apertures;for(const auto& patch:eb.patches){const auto owner=hierarchy.locate_finest(patch.centroid);if(owner.found()&&owner.level==level.level)++patches;}for(int cell:eb.irregular_cells)if(level.owned_cell[cell]){const auto q=eb.grid.cell_coord(cell);const auto box=eb.grid.cell_box(q[0],q[1],q[2]);eb_boxes_.push_back({(float)box.lo.x,(float)box.lo.y,(float)box.lo.z,(float)box.hi.x,(float)box.hi.y,(float)box.hi.z});}for(const auto& problem:eb.unresolved)if(level.owned_cell[problem.parent_cell])++unresolved;}
		const Vec3d size=hierarchy.domain().hi-hierarchy.domain().lo;SimInfo info;info.h=hierarchy.levels().front().h;info.coarse_h=info.h;info.finest_h=hierarchy.finest_cell_size();info.nx=(int)std::llround(size.x/info.h);info.ny=(int)std::llround(size.y/info.h);info.nz=(int)std::llround(size.z/info.h);info.Lx=size.x;info.Ly=size.y;info.Lz=size.z;info.U=run_config.freestream.speed;info.rho=run_config.freestream.rho;info.nu=run_config.freestream.nu;info.name=thin_debug?"thin-y-debug":"paraglider";viewer_->setInfo(info);
		if(thin_debug){viewer_->setMesh(wing);viewer_->setMeshPlacement(ModelPlacement{});thin_debug_display_=true;slice_axis_->setCurrentIndex(1);slice_position_->setValue(500);arrow_mode_->setCurrentIndex(1);tracer_mode_->setCurrentIndex(1);viewer_->setThinDebugState(true,debug_layers);viewer_->setSimulationStep(0);viewer_->frameThinYDebugView();}else{viewer_->setMeshPlacement(config_.placement);viewer_->setThinDebugState(false,0);viewer_->frameWingView();}viewer_->setReferenceU(info.U);updateDebugBoxes();
		const std::size_t bricks=hierarchy.active_brick_count(),dofs=core->pressure_system().storage_size,gpu=core->gpu_bytes();const int muscl=core->embedded_high_order_stencil_count();const QString debug_header=thin_debug?QString("CROPPED-Y DEBUG — %1 CELLS / %2 m\n").arg(debug_layers).arg(debug_width,0,'g',5):QString{};solver_readout_->setText(QString("%1PREPROCESS READY — STARTING CFD\n%2 bricks / %3 pressure DOFs\n%4 EB fragments / %5 apertures / %6 patches\n%7 complete compact MUSCL stencils\n%8 unresolved / %9 static pockets\n%10 sub-grid apertures discarded (%11 m²)\nGPU estimate %12 MiB").arg(debug_header).arg(bricks).arg(dofs).arg(fragments).arg(apertures).arg(patches).arg(muscl).arg(unresolved).arg(static_pockets).arg(discarded_apertures).arg(discarded_aperture_area,0,'g',4).arg(gpu/(1024.0*1024.0),0,'f',1));std::fprintf(stderr,"[paraglider] grid: %zu bricks, %zu DOFs, fragments=%zu apertures=%zu MUSCL=%d patches=%zu unresolved=%zu discarded_apertures=%zu area=%.6g m2, GPU %.2f MiB\n",bricks,dofs,fragments,apertures,muscl,patches,unresolved,discarded_apertures,discarded_aperture_area,gpu/(1024.0*1024.0));spawnWorker(std::move(core));play_button_->setEnabled(true);step_button_->setEnabled(true);if(!play_button_->isChecked())play_button_->setChecked(true);else worker_->setPlaying(true);statusBar()->showMessage(thin_debug?QString("Cropped-Y diagnostic running: %1 cells / %2 m. Watch the STEP counter.").arg(debug_layers).arg(debug_width,0,'g',5):"CFD grid built; paraglider simulation running with freestream +X.",8000);return true;
	}

	void ParagliderWindow::spawnWorker(std::unique_ptr<ExternalAeroCore> core)
	{
		shutdownWorker();snapshot_generation_=surface_generation_=0;worker_=new ParagliderSimWorker(std::move(core));viewer_->setParagliderWorker(worker_);worker_thread_=new QThread(this);worker_->moveToThread(worker_thread_);connect(worker_thread_,&QThread::started,worker_,&ParagliderSimWorker::run);connect(worker_,&ParagliderSimWorker::finished,worker_thread_,&QThread::quit);worker_thread_->start();
	}

	void ParagliderWindow::shutdownWorker()
	{
		if(viewer_)viewer_->setParagliderWorker(nullptr);if(worker_)worker_->stop();if(worker_thread_){worker_thread_->quit();worker_thread_->wait();}delete worker_;worker_=nullptr;delete worker_thread_;worker_thread_=nullptr;snapshot_generation_=surface_generation_=0;
	}

	void ParagliderWindow::applySurfaceColour()
	{
		if(!viewer_||delta_cp_.empty())return;const int mode=surface_colour_->currentIndex();if(mode==0)viewer_->setTriangleSurfaceCp(cp_plus_,cp_minus_,-side_cp_range_,side_cp_range_);else if(mode==2)viewer_->setTriangleSurfaceCp(cp_plus_,cp_plus_,-side_cp_range_,side_cp_range_);else if(mode==3)viewer_->setTriangleSurfaceCp(cp_minus_,cp_minus_,-side_cp_range_,side_cp_range_);else viewer_->setTriangleSurfaceCp(delta_cp_,delta_cp_,-delta_cp_range_,delta_cp_range_);
	}

	void ParagliderWindow::updateDebugBoxes(){if(viewer_)viewer_->setParagliderDebugBoxes(show_amr_&&show_amr_->isChecked()?amr_boxes_:std::vector<std::array<float,6>>{},show_eb_&&show_eb_->isChecked()?eb_boxes_:std::vector<std::array<float,6>>{});}

	void ParagliderWindow::updateSnapshot()
	{
		if(!worker_)return;ParagliderDisplaySnapshot s;if(!worker_->latestSnapshot(snapshot_generation_,surface_generation_,s))return;if(!s.error.empty()){solver_readout_->setText(QString("SOLVER STOPPED: %1").arg(QString::fromStdString(s.error)));play_button_->setChecked(false);return;}if(!s.initialized){solver_readout_->setText("Initializing pressure field on GPU...");return;}last_steps_=s.steps;last_time_=s.physical_time;viewer_->setSimulationStep(s.steps);if(!s.delta_cp.empty()){cp_plus_=std::move(s.cp_plus);cp_minus_=std::move(s.cp_minus);delta_cp_=std::move(s.delta_cp);delta_cp_range_=std::max(std::abs(s.cp_min),std::abs(s.cp_max));side_cp_range_=std::max(std::abs(s.side_cp_min),std::abs(s.side_cp_max));applySurfaceColour();}solver_readout_->setText(QString("step %1   t=%2 s   dt=%3 ms   CFL=%4\nGPU step %5 ms; projection %6 ms / %7 it\nresidual %8; regular/EB max %9 / %10 m/s\nEB transport rate %11 1/s; GPU %12 MiB").arg(s.steps).arg(s.physical_time,0,'f',4).arg(s.dt*1e3,0,'f',3).arg(s.effective_cfl,0,'f',2).arg(s.step_ms,0,'f',1).arg(s.projection_ms,0,'f',1).arg(s.pressure_iterations).arg(s.residual,0,'g',3).arg(s.max_abs_regular_velocity,0,'g',4).arg(s.max_abs_special_velocity,0,'g',4).arg(s.max_embedded_cfl_rate,0,'g',4).arg(s.gpu_bytes/(1024.0*1024.0),0,'f',1));QString loads=QString("PRESSURE-ONLY LOADS\nD +X %1 N   S +Y %2 N   L +Z %3 N\nΔCp %4 ... %5").arg(s.pressure_force.x,0,'f',3).arg(s.pressure_force.y,0,'f',3).arg(s.pressure_force.z,0,'f',3).arg(s.cp_min,0,'f',3).arg(s.cp_max,0,'f',3);if(s.coefficients_valid)loads+=QString("\nCd,p %1   Cs,p %2   Cl,p %3").arg(s.cd_pressure,0,'f',4).arg(s.cs_pressure,0,'f',4).arg(s.cl_pressure,0,'f',4);else loads+="\nCL/CD withheld: reference area is zero";if(s.conservation_valid)loads+=QString("\nflux error max/sum/net %1 / %2 / %3 m³/s").arg(s.max_integrated_flux_error,0,'g',3).arg(s.absolute_integrated_flux_error,0,'g',3).arg(s.net_integrated_flux_error,0,'g',3);if(s.max_abs_special_velocity>std::max(50.0,3*s.max_abs_regular_velocity))loads+="\nWARNING: compact-EB hotspot; loads not trustworthy";loads+="\nskin friction not implemented";load_readout_->setText(loads);
	}

	void ParagliderWindow::closeEvent(QCloseEvent* event){shutdownWorker();QMainWindow::closeEvent(event);}
}
