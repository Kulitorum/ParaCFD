#include "core/cuda_probe.h"
#include "gui/paraglider_window.h"
#include "gui/slice_viewer.h"

#include <QApplication>
#include <QImage>
#include <QSurfaceFormat>
#include <QTimer>

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>

int main(int argc,char** argv)
{
	std::string config,step,screenshot,window_screenshot,visualization_preset;int autoclose_ms=0;double thin_y_fraction=-1,thin_y_width=0.125,sweep_min=0,sweep_max=0,sweep_step=1,sweep_mean_tolerance_percent=2.0,sweep_max_flow_throughs=2.5,reference_area=-1;bool offscreen=false,thin_y_width_explicit=false,conservative_momentum=false,half_wing=false,sweep_requested=false,exit_after_sweep=false,exit_after_auto_pause=false;

	for(int i=1;i<argc;++i)
	{
		const std::string argument=argv[i];
		if(argument=="--config"&&i+1<argc)config=argv[++i];
		else if(argument=="--load-step"&&i+1<argc)step=argv[++i];
		else if(argument=="--autoclose-ms"&&i+1<argc)autoclose_ms=std::atoi(argv[++i]);
		else if(argument=="--screenshot"&&i+1<argc)screenshot=argv[++i];
		else if(argument=="--window-screenshot"&&i+1<argc)window_screenshot=argv[++i];
		else if(argument=="--visualization-preset"&&i+1<argc)visualization_preset=argv[++i];
		else if(argument=="--thin-y-fraction"&&i+1<argc)thin_y_fraction=std::atof(argv[++i]);
		else if(argument=="--thin-y-width"&&i+1<argc){thin_y_width=std::atof(argv[++i]);thin_y_width_explicit=true;}
		else if(argument=="--conservative-momentum")conservative_momentum=true;
		else if(argument=="--staggered-momentum")conservative_momentum=false;
		else if(argument=="--half-wing")half_wing=true;
		else if(argument=="--aoa-sweep"&&i+3<argc){sweep_min=std::atof(argv[++i]);sweep_max=std::atof(argv[++i]);sweep_step=std::atof(argv[++i]);sweep_requested=true;}
		else if(argument=="--sweep-mean-tolerance-percent"&&i+1<argc)sweep_mean_tolerance_percent=std::atof(argv[++i]);
		else if(argument=="--sweep-max-flow-throughs"&&i+1<argc)sweep_max_flow_throughs=std::atof(argv[++i]);
		else if(argument=="--reference-area"&&i+1<argc)reference_area=std::atof(argv[++i]);
		else if(argument=="--exit-after-sweep")exit_after_sweep=true;
		else if(argument=="--exit-after-auto-pause")exit_after_auto_pause=true;
		else if(argument=="--offscreen")offscreen=true;
		else if(!argument.empty()&&argument[0]!='-')config=argument;
	}
	if(offscreen)qputenv("QT_QPA_PLATFORM","offscreen");
	QSurfaceFormat format;format.setDepthBufferSize(24);format.setSamples(4);format.setProfile(QSurfaceFormat::CoreProfile);format.setVersion(4,3);QSurfaceFormat::setDefaultFormat(format);
	QApplication app(argc,argv);QApplication::setApplicationName("ParaCFD");QApplication::setOrganizationName("ParaCFD");
	std::string error,gpu;int major=0,minor=0;if(paracfd::core::cuda_probe_usable(&gpu,&major,&minor,&error)!=0){std::fprintf(stderr,"[paraglider] CUDA unavailable: %s\n",error.c_str());return 2;}std::fprintf(stderr,"[paraglider] CUDA device: %s (compute %d.%d)\n",gpu.c_str(),major,minor);
	paracfd::gui::ParagliderWindow window;window.show();if(!visualization_preset.empty())window.setVisualizationPreset(QString::fromStdString(visualization_preset));if(thin_y_fraction>=0||thin_y_width_explicit)window.setThinYDiagnostic(thin_y_fraction>=0?thin_y_fraction:0.5,thin_y_width);window.setConservativeMomentum(conservative_momentum);window.setAoaSweepControls(sweep_mean_tolerance_percent,sweep_max_flow_throughs,reference_area);bool ready=false;if(!config.empty()){ready=window.loadConfigFile(QString::fromStdString(config),false);if(half_wing)window.setHalfWingSimulation(true);window.setAoaSweepControls(sweep_mean_tolerance_percent,sweep_max_flow_throughs,reference_area);if(ready)ready=sweep_requested?window.startAoaSweep(sweep_min,sweep_max,sweep_step):window.buildGrid();}else if(!step.empty()){if(half_wing)window.setHalfWingSimulation(true);ready=window.loadStepFile(QString::fromStdString(step),true);window.setAoaSweepControls(sweep_mean_tolerance_percent,sweep_max_flow_throughs,reference_area);if(ready)ready=sweep_requested?window.startAoaSweep(sweep_min,sweep_max,sweep_step):window.buildGrid();}else{window.loadLastFile();if(half_wing)window.setHalfWingSimulation(true);window.setAoaSweepControls(sweep_mean_tolerance_percent,sweep_max_flow_throughs,reference_area);if(sweep_requested)ready=window.startAoaSweep(sweep_min,sweep_max,sweep_step);}
	QTimer sweep_exit_timer;if(exit_after_sweep&&sweep_requested&&ready){QObject::connect(&sweep_exit_timer,&QTimer::timeout,&app,[&]{if(!window.aoaSweepActive())app.quit();});sweep_exit_timer.start(250);}else if(exit_after_sweep&&sweep_requested)QTimer::singleShot(0,&app,&QApplication::quit);
	bool auto_pause_exit_failed=false;QTimer auto_pause_exit_timer;
	if(exit_after_auto_pause&&!sweep_requested&&ready)
	{
		QObject::connect(&auto_pause_exit_timer,&QTimer::timeout,&app,[&]
		{
			if(window.simulationAutoPaused())app.quit();
			else if(!window.simulationRunning()){auto_pause_exit_failed=true;app.quit();}
		});
		auto_pause_exit_timer.start(250);
	}
	else if(exit_after_auto_pause&&!sweep_requested)QTimer::singleShot(0,&app,&QApplication::quit);
	if(autoclose_ms>0){const int capture_delay=std::max(100,autoclose_ms-500);if(!screenshot.empty())QTimer::singleShot(capture_delay,&app,[&window,screenshot]{if(window.viewer())window.viewer()->grabFramebuffer().save(QString::fromStdString(screenshot));});if(!window_screenshot.empty())QTimer::singleShot(capture_delay,&app,[&window,window_screenshot]{window.grab().save(QString::fromStdString(window_screenshot));});QTimer::singleShot(autoclose_ms,&app,&QApplication::quit);}
	const int result=app.exec();const bool expected_auto_pause_missing=exit_after_auto_pause&&!sweep_requested&&ready&&!window.simulationAutoPaused();window.close();const double seconds=autoclose_ms/1000.0;std::fprintf(stderr,"[paraglider] exit: %lld steps, t=%.6f s, %lld frames%s\n",window.steps(),window.physicalTime(),window.viewer()?window.viewer()->framesRendered():0,autoclose_ms>0?" (scripted smoke)":"");if((autoclose_ms>0&&ready&&window.steps()==0)||(exit_after_sweep&&sweep_requested&&!ready)||(exit_after_auto_pause&&!sweep_requested&&!ready)||auto_pause_exit_failed||expected_auto_pause_missing)return 3;(void)seconds;return result;

}
