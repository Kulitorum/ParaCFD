// gl_thread_check.h — G1 hard rule: every OpenGL call happens on the main (GUI) thread.
// PARACFD_ASSERT_GL_THREAD() fires a Qt debug assertion (and aborts in a debug build) if a
// GL entry point is ever reached off the main thread. This is the automated half of the
// PLAN §4 "no GL calls off the main thread (assert via Qt debug)" gate.
#pragma once

#include <QCoreApplication>
#include <QThread>
#include <QtGlobal>

namespace paracfd::gui
{
	inline bool on_main_thread()
	{
		QCoreApplication* app = QCoreApplication::instance();
		return app && QThread::currentThread() == app->thread();
	}
}

#define PARACFD_ASSERT_GL_THREAD() \
	Q_ASSERT_X(paracfd::gui::on_main_thread(), Q_FUNC_INFO, "OpenGL call issued off the main thread")
