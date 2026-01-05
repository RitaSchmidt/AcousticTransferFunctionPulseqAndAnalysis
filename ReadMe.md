# Measure Acoustic Transfer Functions

# An open sequence (via Pulseq) and analysis for measuring acoustic/vibrational transfer functions of MRI systems

## Introduction

During an MRI scan, the gradient coils (responsible for encoding position) are driven by currents through them. These currents run within the large static magnetic field of the MRI, B0, which exerts a Lorentz force on the wires of the gradient coils. As the currents change during the scan, the Lorentz force also changes, leading to vibrations of the gradient coils. (Other sources of vibrations also exist; any induced eddy current in conducting parts will also result in Lorentz forces and vibrations.) These vibrations lead to the typical sounds of the MRI scans.

The mechanical response of the system to these gradient induced vibrations depends on their component frequencies. Frequencies close to mechanical resonances will be amplified relative to others. Such frequencies make the scans loud, requiring the subjects to wear protective gear (earplugs and/or earphones), and in some cases the vibrations are so strong that they may damage the system itself.

Therefore, it is desirable to know the response of MRI systems to different frequencies. This response is known as the transfer function. The code here should allow to perform such measurements safely, without damaging the system. However, ***no warranty is given and you should know what you are doing.***

**NOTE:** From our experience, transfer function tend to change over time (and even within a session).

## General Usage

The script `GenerateAcousticTrasnferFunctionScan.m` produces a Pulseq `.seq` file driving one of the gradients (X, Y, or Z) with a series of equispaced (modulated) triangles. See the script's help and `Inputs` section for the different parameters and their meaning.

While the scan is running the resulting sounds or vibrations should be recorded (with MRI compatible equipment) and saved. Supported formats are any format Matlab's `audioread()` supports (`.aifc`, `.aiff`, `.aif`, `.au`, `.flac`, `.ogg`, `.opus`, `.wav`, `.mp3`, `.m4a`, `.mp4`), as well as CSV and Matlab's time tables inside `.mat` files. CSV files can be either comma separated (as the name suggest), or semicolon separated (in case the comma is used to mark the decimal separator). Note that a header of length 16 lines is expected, but you can edit `CSVNumHeaderLinesDefault` if your header differs. In case of a Matlab `.mat` file, it is expected to include only a single timetable variable and that table should include two columns, one of which is `Time` (case insensitive).

Once the recordings are done, update the paths in the `RunRecording2TransferFunction.m` script, which will pass the recording and the Pulseq `.seq` file to the `Recording2TransferFunction.m` which will generate the transfer function (a `.mat` file).

The resulting `.mat` file will include a structure three parameters defining the transfer function (names are kept for backwards compatibility, but you may change them):
* `FreqsRaster` - The frequencies (Hz) in which the transfer function is sampled.
* `RRaster` - The transfer function itself (at the frequency samples `FreqsRaster`). The units are arbitrary and different transfer functions, in general, cannot be compared, except for the relative change within them.
* `RErrRaster` - The standard deviation of each transfer function sample. It is based on the averaging done during the measurement itself.

### Sample results

Here is a zoom between 300 Hz and 1500 Hz, of a transfer function we measured.

![Sample, zoomed in, tranfer function](https://github.com/RitaSchmidt/AcousticTransferFunctionPulseqAndAnalysis/blob/main/SampleData/AnalysisOutput/Amplification_vs_Freq.svg?raw=true)
