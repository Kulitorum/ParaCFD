// video_recorder.cpp — see video_recorder.h. Streams raw RGBA frames to an ffmpeg subprocess.
#include "gui/video_recorder.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QProcess>
#include <QStandardPaths>
#include <QStringList>

#include <cstdio>

namespace windcfd::gui
{
	VideoRecorder::VideoRecorder() = default;
	VideoRecorder::~VideoRecorder() { finish(); }

	QString VideoRecorder::findFfmpeg()
	{
		// 1) explicit override.
		const QByteArray env = qgetenv("WINDCFD_FFMPEG");
		if (!env.isEmpty())
		{
			const QString p = QString::fromLocal8Bit(env);
			if (QFileInfo::exists(p)) return p;
		}
		// 2) bundled next to the executable — the installer ships ffmpeg.exe into {app}, so a deployed
		// install records video out of the box without ffmpeg on PATH.
		if (QCoreApplication::instance())
		{
			const QString bundled = QDir(QCoreApplication::applicationDirPath()).filePath("ffmpeg.exe");
			if (QFileInfo::exists(bundled)) return bundled;
		}
		// 3) on PATH.
		const QString onpath = QStandardPaths::findExecutable("ffmpeg");
		if (!onpath.isEmpty()) return onpath;
		// 4) the self-contained build present on this dev machine (last-resort dev convenience).
		const QString known = "C:/Users/Micro/Downloads/ffmpeg-2025-12-22-git-c50e5c7778-essentials_build/bin/ffmpeg.exe";
		if (QFileInfo::exists(known)) return known;
		return QString();
	}

	bool VideoRecorder::start(const QString& path, int crf, const QString& preset, int fps)
	{
		if (recording_) return true;
		ffmpeg_ = findFfmpeg();
		if (ffmpeg_.isEmpty())
		{
			std::fprintf(stderr, "[video] ffmpeg not found (set WINDCFD_FFMPEG, or put ffmpeg on PATH); recording disabled\n");
			return false;
		}
		path_ = path;
		crf_ = (crf < 0) ? 0 : (crf > 51) ? 51 : crf;                 // valid x264 CRF band
		fps_ = (fps < 1) ? 1 : fps;
		preset_ = preset.isEmpty() ? QStringLiteral("fast") : preset;
		frames_written_ = 0;
		recording_ = true; // armed; ffmpeg spawns on the first frame (locks the size)
		started_ = false;
		std::fprintf(stderr, "[video] recording armed -> %s (crf %d, preset %s, %d fps)\n",
			path_.toUtf8().constData(), crf_, preset_.toUtf8().constData(), fps_);
		return true;
	}

	void VideoRecorder::writeFrame(const QImage& frame)
	{
		if (!recording_ || frame.isNull()) return;

		// Spawn ffmpeg on the first frame, locking the output size to this frame (even dims for yuv420p).
		if (!started_)
		{
			w_ = frame.width() & ~1;
			h_ = frame.height() & ~1;
			if (w_ < 2 || h_ < 2) return; // wait for a real-sized frame
			proc_ = std::make_unique<QProcess>();
			proc_->setProcessChannelMode(QProcess::ForwardedErrorChannel); // ffmpeg's log -> our stderr
			const QStringList args = {
				"-y",
				"-f", "rawvideo",
				"-pixel_format", "rgba",
				"-video_size", QStringLiteral("%1x%2").arg(w_).arg(h_),
				"-framerate", QString::number(fps_),
				"-i", "pipe:0",
				"-an",
				"-c:v", "libx264",
				"-crf", QString::number(crf_),
				"-preset", preset_, // crf sets the quality; the preset trades encode speed vs size
				"-pix_fmt", "yuv420p",
				"-movflags", "+faststart",
				path_
			};
			proc_->start(ffmpeg_, args);
			if (!proc_->waitForStarted(5000))
			{
				std::fprintf(stderr, "[video] ffmpeg failed to start; recording disabled\n");
				proc_.reset();
				recording_ = false;
				return;
			}
			started_ = true;
			std::fprintf(stderr, "[video] ffmpeg streaming %dx%d (h264 crf%d preset %s @ %d fps) -> %s\n",
				w_, h_, crf_, preset_.toUtf8().constData(), fps_, path_.toUtf8().constData());
		}

		if (!proc_ || proc_->state() != QProcess::Running) return;

		// Backpressure: if ffmpeg is more than ~3 frames behind, drop this one rather than stall the UI.
		if (proc_->bytesToWrite() > (qint64)w_ * h_ * 4 * 3) return;

		// Tightly-packed RGBA rows (Format_RGBA8888 is 4 B/px ⇒ no per-row padding), scaled to the lock.
		QImage img = frame.convertToFormat(QImage::Format_RGBA8888);
		if (img.width() != w_ || img.height() != h_)
			img = img.scaled(w_, h_, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);

		const qsizetype rowbytes = (qsizetype)w_ * 4;
		if (img.bytesPerLine() == rowbytes)
			proc_->write(reinterpret_cast<const char*>(img.constBits()), rowbytes * h_);
		else
			for (int y = 0; y < h_; ++y)
				proc_->write(reinterpret_cast<const char*>(img.constScanLine(y)), rowbytes);
		++frames_written_;
	}

	void VideoRecorder::flush()
	{
		if (!proc_) return;
		while (proc_->bytesToWrite() > 0)
			if (!proc_->waitForBytesWritten(3000)) break;
	}

	void VideoRecorder::finish()
	{
		if (proc_)
		{
			proc_->closeWriteChannel(); // EOF on stdin -> ffmpeg flushes + finalizes the MP4
			if (!proc_->waitForFinished(20000))
			{
				proc_->kill();
				proc_->waitForFinished(2000);
				std::fprintf(stderr, "[video] ffmpeg did not finish in time; killed (file may be truncated)\n");
			}
			else
				std::fprintf(stderr, "[video] finished -> %s (ffmpeg exit %d)\n", path_.toUtf8().constData(), proc_->exitCode());
			proc_.reset();
		}
		recording_ = false;
		started_ = false;
	}
}
