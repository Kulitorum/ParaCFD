// video_settings_dialog.h — non-modal settings popup for the MP4 recorder (File → "Record Video…").
//
// Exposes the recording knobs that were previously hard-coded or CLI-only: output path, capture cadence
// (sim-steps per frame), H.264 quality (CRF + x264 preset) and playback fps. The dialog owns NO recording
// logic — it emits startRequested()/stopRequested() and MainWindow reads the getters, starts/stops the
// VideoRecorder, and pushes the live frame count back via setRecordingStatus(). Non-modal so the sim keeps
// running while it is open.
#pragma once

#include <QDialog>

class QComboBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;

namespace scour::gui
{
	class VideoSettingsDialog : public QDialog
	{
		Q_OBJECT
	public:
		explicit VideoSettingsDialog(QWidget* parent = nullptr);

		// Seed the controls from the caller's current values (call before show()).
		void setValues(const QString& path, int cadenceSteps, int crf, const QString& preset, int fps);

		// Chosen values (read by MainWindow on startRequested()).
		QString path() const;
		int cadenceSteps() const;
		int crf() const;
		QString preset() const;
		int fps() const;

		// Reflect the live recorder state: flips the buttons + shows "● REC — N frames" / "idle".
		void setRecordingStatus(bool recording, long long framesWritten);

	signals:
		void startRequested();
		void stopRequested();

	private:
		void browse();

		QLineEdit* path_edit_ = nullptr;
		QPushButton* browse_btn_ = nullptr;
		QSpinBox* cadence_spin_ = nullptr;
		QSpinBox* crf_spin_ = nullptr;
		QComboBox* preset_box_ = nullptr;
		QSpinBox* fps_spin_ = nullptr;
		QPushButton* start_btn_ = nullptr;
		QPushButton* stop_btn_ = nullptr;
		QLabel* status_ = nullptr;
	};
}
