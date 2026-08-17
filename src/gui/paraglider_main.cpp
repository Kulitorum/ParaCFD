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
	std::string config,step,screenshot;int autoclose_ms=0;double thin_y_fraction=-1,thin_y_width=0.125;bool offscreen=false,thin_y_width_explicit=false,conservative_momentum=true;
	for(int i=1;i<argc;++i)
	{
		const std::string argument=argv[i];
		if(argument=="--config"&&i+1<argc)config=argv[++i];
		else if(argument=="--load-step"&&i+1<argc)step=argv[++i];
		else if(argument=="--autoclose-ms"&&i+1<argc)autoclose_ms=std::atoi(argv[++i]);
		else if(argument=="--screenshot"&&i+1<argc)screenshot=argv[++i];
		else if(argument=="--thin-y-fraction"&&i+1<argc)thin_y_fraction=std::atof(argv[++i]);
		else if(argument=="--thin-y-width"&&i+1<argc){thin_y_width=std::atof(argv[++i]);thin_y_width_explicit=true;}
		else if(argument=="--conservative-momentum")conservative_momentum=true;
		else if(argument=="--staggered-momentum")conservative_momentum=false;
		else if(argument=="--offscreen")offscreen=true;
		else if(!argument.empty()&&argument[0]!='-')config=argument;
	}
	if(offscreen)qputenv("QT_QPA_PLATFORM","offscreen");
	QSurfaceFormat format;format.setDepthBufferSize(24);format.setSamples(4);format.setProfile(QSurfaceFormat::CoreProfile);format.setVersion(4,3);QSurfaceFormat::setDefaultFormat(format);
	QApplication app(argc,argv);QApplication::setApplicationName("ParaCFD");QApplication::setOrganizationName("ParaCFD");
	std::string error,gpu;int major=0,minor=0;if(paracfd::core::cuda_probe_usable(&gpu,&major,&minor,&error)!=0){std::fprintf(stderr,"[paraglider] CUDA unavailable: %s\n",error.c_str());return 2;}std::fprintf(stderr,"[paraglider] CUDA device: %s (compute %d.%d)\n",gpu.c_str(),major,minor);
	paracfd::gui::ParagliderWindow window;window.show();if(thin_y_fraction>=0||thin_y_width_explicit)window.setThinYDiagnostic(thin_y_fraction>=0?thin_y_fraction:0.5,thin_y_width);window.setConservativeMomentum(conservative_momentum);bool ready=false;if(!config.empty())ready=window.loadConfigFile(QString::fromStdString(config),true);else if(!step.empty()){ready=window.loadStepFile(QString::fromStdString(step),true)&&window.buildGrid();}else window.loadLastFile();
	if(autoclose_ms>0){if(!screenshot.empty())QTimer::singleShot(std::max(100,autoclose_ms-500),&app,[&window,screenshot]{if(window.viewer())window.viewer()->grabFramebuffer().save(QString::fromStdString(screenshot));});QTimer::singleShot(autoclose_ms,&app,&QApplication::quit);}
	const int result=app.exec();window.close();const double seconds=autoclose_ms/1000.0;std::fprintf(stderr,"[paraglider] exit: %lld steps, t=%.6f s, %lld frames%s\n",window.steps(),window.physicalTime(),window.viewer()?window.viewer()->framesRendered():0,autoclose_ms>0?" (scripted smoke)":"");if(autoclose_ms>0&&ready&&window.steps()==0)return 3;(void)seconds;return result;
}
