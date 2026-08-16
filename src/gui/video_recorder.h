// video_recorder.h — pipe rendered GL frames to ffmpeg for a high-quality MP4 of the run.
//
// paracfd-gui captures the viewer's framebuffer whenever the sand bed changes and streams the frames
// as raw RGBA into an ffmpeg subprocess (H.264 / crf 17). No libav linkage — ffmpeg.exe is spawned
// via QProcess, so there is nothing to build/deploy and the codec licence stays with the exe.
//
// Everything here runs on the MAIN (GL) thread: the caller grabs the framebuffer (QOpenGLWidget::
// grabFramebuffer) and hands each QImage to writeFrame(). Writes are non-blocking (Qt buffers them,
// the event loop flushes to the pipe); a frame is dropped if ffmpeg falls behind, so the UI never
// stalls. finish() closes ffmpeg's stdin and waits for the MP4 to be muxed + closed.
#pragma once

#include <QImage>
#include <QString>

#include <memory>

class QProcess;

namespace paracfd::gui
{
	class VideoRecorder
	{
	public:
		VideoRecorder();
		~VideoRecorder(); // finalizes the file (idempotent) so the MP4 is always closed on destruction

		VideoRecorder(const VideoRecorder&) = delete;
		VideoRecorder& operator=(const VideoRecorder&) = delete;

		// Arm recording to `path` (an .mp4). ffmpeg is located (PARACFD_FFMPEG env → PATH → bundled build)
		// and spawned LAZILY on the first frame, which locks the output size to that frame's dimensions
		// (rounded down to even, as yuv420p requires). `crf` is the H.264 quality (0 lossless … 51 worst,
		// default 17 visually lossless), `preset` the x264 speed/size trade (ultrafast…veryslow, default
		// "fast"), `fps` the playback frame rate the frames are muxed at (default 30). Returns false (and
		// stays disabled) if ffmpeg is not found or is already recording.
		bool start(const QString& path, int crf = 17, const QString& preset = QStringLiteral("fast"), int fps = 30);

		// Write one frame. Any size/format is accepted: it is converted to RGBA8888 and scaled to the
		// locked size. No-op when not recording; drops the frame (no stall) if the pipe is backed up.
		void writeFrame(const QImage& frame);

		// Block until the buffered frame data has been written to ffmpeg. The live GUI never needs this
		// (its event loop drains the pipe between captures); it exists for the headless self-test, which
		// has no event loop, so the non-blocking writes would otherwise back up and be dropped.
		void flush();

		// Close ffmpeg's stdin and wait for it to finish writing the MP4. Idempotent; safe on close.
		void finish();

		bool recording() const { return recording_; }
		QString path() const { return path_; }
		// Number of frames actually piped to ffmpeg since the last start() (dropped/backpressured frames are
		// NOT counted). Drives the live "● REC — N frames" status in the recording dialog.
		long long framesWritten() const { return frames_written_; }

	private:
		static QString findFfmpeg();

		std::unique_ptr<QProcess> proc_;
		QString ffmpeg_;         // resolved ffmpeg executable
		QString path_;           // output .mp4
		bool recording_ = false; // armed or streaming
		bool started_ = false;   // ffmpeg spawned (size locked)
		int w_ = 0, h_ = 0;      // locked, even frame size
		int crf_ = 17;           // H.264 quality (0 lossless … 51 worst); set by start()
		int fps_ = 30;           // playback frame rate the frames are muxed at; set by start()
		QString preset_ = QStringLiteral("fast"); // x264 speed/size preset; set by start()
		long long frames_written_ = 0;            // frames piped to ffmpeg since start()
	};
}
