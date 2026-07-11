// video_settings_dialog.cpp — see video_settings_dialog.h.
#include "gui/video_settings_dialog.h"

#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

namespace windcfd::gui
{
	// x264 speed/size presets, slowest→fastest encode is the reverse; slower = smaller file at the same CRF.
	static const char* const kPresets[] = {
		"ultrafast", "superfast", "veryfast", "faster", "fast", "medium", "slow", "slower", "veryslow"
	};

	VideoSettingsDialog::VideoSettingsDialog(QWidget* parent)
		: QDialog(parent)
	{
		setWindowTitle("Record Video — settings");
		setModal(false); // non-modal: the sim keeps running while this is open

		QVBoxLayout* col = new QVBoxLayout(this);
		QFormLayout* form = new QFormLayout;
		form->setLabelAlignment(Qt::AlignLeft);

		// Output file + Browse.
		QHBoxLayout* pathRow = new QHBoxLayout;
		path_edit_ = new QLineEdit;
		path_edit_->setPlaceholderText("scour_run.mp4");
		path_edit_->setToolTip("Output .mp4 path.");
		browse_btn_ = new QPushButton("Browse…");
		connect(browse_btn_, &QPushButton::clicked, this, [this] { browse(); });
		pathRow->addWidget(path_edit_, 1);
		pathRow->addWidget(browse_btn_);
		form->addRow("File", pathRow);

		// Capture cadence — one frame is captured every N simulation steps. It is a TIME-LAPSE: lower = more
		// frames = smoother but a longer clip; higher = fewer frames = a shorter, faster-moving clip.
		cadence_spin_ = new QSpinBox;
		cadence_spin_->setRange(1, 100000);
		cadence_spin_->setSingleStep(5);
		cadence_spin_->setSuffix(" steps/frame");
		cadence_spin_->setToolTip("Capture one frame every N simulation steps (a time-lapse). Lower = more, smoother frames but a longer clip; higher = fewer frames, a shorter clip.");
		form->addRow("Cadence", cadence_spin_);

		// CRF — H.264 quality: 0 lossless … 51 worst; 17 ≈ visually lossless. Lower = better + bigger.
		crf_spin_ = new QSpinBox;
		crf_spin_->setRange(0, 51);
		crf_spin_->setToolTip("H.264 quality (CRF). 0 = lossless, 17 ≈ visually lossless, 23 = default, 51 = worst. Lower = better quality and a bigger file.");
		form->addRow("Quality (CRF)", crf_spin_);

		// Preset — x264 speed/size trade at a fixed CRF. Slower = smaller file, more CPU.
		preset_box_ = new QComboBox;
		for (const char* p : kPresets) preset_box_->addItem(p);
		preset_box_->setToolTip("x264 encode preset. Slower presets shrink the file at the same quality but use more CPU; 'fast' keeps encoding off the render thread's back.");
		form->addRow("Preset", preset_box_);

		// Playback fps — the rate the captured frames are muxed at (independent of the capture cadence).
		fps_spin_ = new QSpinBox;
		fps_spin_->setRange(1, 120);
		fps_spin_->setSuffix(" fps");
		fps_spin_->setToolTip("Playback frame rate the captured frames are muxed at (independent of the capture cadence).");
		form->addRow("Playback fps", fps_spin_);

		col->addLayout(form);

		QFrame* line = new QFrame; line->setFrameShape(QFrame::HLine); line->setEnabled(false);
		col->addWidget(line);

		// Live status + Start/Stop.
		status_ = new QLabel("idle");
		status_->setToolTip("Recording state and the number of frames written so far.");
		col->addWidget(status_);
		QHBoxLayout* btnRow = new QHBoxLayout;
		start_btn_ = new QPushButton("Start");
		stop_btn_ = new QPushButton("Stop");
		connect(start_btn_, &QPushButton::clicked, this, [this] { emit startRequested(); });
		connect(stop_btn_, &QPushButton::clicked, this, [this] { emit stopRequested(); });
		btnRow->addWidget(start_btn_);
		btnRow->addWidget(stop_btn_);
		btnRow->addStretch(1);
		QPushButton* closeBtn = new QPushButton("Close");
		connect(closeBtn, &QPushButton::clicked, this, [this] { hide(); });
		btnRow->addWidget(closeBtn);
		col->addLayout(btnRow);

		setRecordingStatus(false, 0); // initial button/label state
	}

	void VideoSettingsDialog::browse()
	{
		const QString start = path_edit_->text().isEmpty() ? QStringLiteral("scour_run.mp4") : path_edit_->text();
		const QString fn = QFileDialog::getSaveFileName(this, "Record Video", start, "MP4 video (*.mp4)");
		if (!fn.isEmpty()) path_edit_->setText(fn);
	}

	void VideoSettingsDialog::setValues(const QString& path, int cadenceSteps, int crf, const QString& preset, int fps)
	{
		if (!path.isEmpty()) path_edit_->setText(path);
		cadence_spin_->setValue(cadenceSteps);
		crf_spin_->setValue(crf);
		int pi = preset_box_->findText(preset);
		preset_box_->setCurrentIndex(pi >= 0 ? pi : preset_box_->findText("fast"));
		fps_spin_->setValue(fps);
	}

	QString VideoSettingsDialog::path() const { return path_edit_->text(); }
	int VideoSettingsDialog::cadenceSteps() const { return cadence_spin_->value(); }
	int VideoSettingsDialog::crf() const { return crf_spin_->value(); }
	QString VideoSettingsDialog::preset() const { return preset_box_->currentText(); }
	int VideoSettingsDialog::fps() const { return fps_spin_->value(); }

	void VideoSettingsDialog::setRecordingStatus(bool recording, long long framesWritten)
	{
		start_btn_->setEnabled(!recording);
		stop_btn_->setEnabled(recording);
		// The settings are locked while streaming (the size/quality is fixed on the first frame).
		path_edit_->setEnabled(!recording);
		browse_btn_->setEnabled(!recording);
		cadence_spin_->setEnabled(!recording);
		crf_spin_->setEnabled(!recording);
		preset_box_->setEnabled(!recording);
		fps_spin_->setEnabled(!recording);
		if (recording)
		{
			status_->setText(QString("<b><span style='color:#e44;'>●</span> REC</b> — %1 frames").arg(framesWritten));
		}
		else
			status_->setText("idle");
	}
}
