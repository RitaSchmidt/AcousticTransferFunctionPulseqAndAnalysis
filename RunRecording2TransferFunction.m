%% Setup

% "Base" folder to look for sub-folders and files
BaseFolder = fullfile(pwd, 'SampleData') ;

% Full path to audio recording file
AudioFilePath = fullfile(BaseFolder, 'Recordings', ...
                         'AudioRecording-2025.12.09-13_13.wav') ;
% Full path to .seq file matching the audio recording.
PulseqSeqFilePath = fullfile(BaseFolder, 'PulseqFiles', ...
                             'ScanAcousticTransferFunc-Reps10.seq') ;


% Set full path to output file (will hold the "transfer function").
OutputFolder = fullfile(BaseFolder, 'AnalysisOutput') ;
[~, AudioFileName, AudioFileExt] = fileparts(AudioFilePath) ;
OutpufFile = fullfile(OutputFolder, ['Amplification-', AudioFileName '.mat']) ;

% Define range of frequencies to display
FreqRangeShow = [300, 1500] ;

% Use logarithmic scaling of transfer function plot?
bPlotLogScale = true ;

% How many initial "TRs" to discard, to allow for acoustic "steady state".
TR2AcousticSteadyState = 1 ; 

% Whether to save the results in a .mat file
bSaveOutput = true ;


%% Determine trasnfer function

% Read audio-file, matching .seq file, and determine transfer function.
[TransferFuncStruct] = Recording2TransferFunction(AudioFilePath, ...
                                                  PulseqSeqFilePath, ...
                                                  TR2AcousticSteadyState) ;

% % Full set of parameters/inputs for Recording2TransferFunction()
% TR2AcousticSteadyState = [] ;
% TransferFuncResHz = [] ; % read from .seq file
% TransferFuncPosBWHz = [] ; % read from .seq file
% AudioFormat = 'CSV' ;
% CSVNumHeaderLines = 16 ;
% bDecimalPointCSV = true ;
% 
% [TransferFuncStruct] = Recording2TransferFunction(AudioFilePath, ...
%                                                   PulseqSeqFilePath, ...
%                                                   TR2AcousticSteadyState, ...
%                                                   TransferFuncResHz, ...
%                                                   TransferFuncPosBWHz, ...
%                                                   AudioFormat, ...
%                                                   CSVNumHeaderLines, ...
%                                                   bDecimalPointCSV) ;

%% Plot results:

figure ;
  errorbar(TransferFuncStruct.f, TransferFuncStruct.Amplification, ...
           TransferFuncStruct.AmplificationStdErr) ;
  xlabel('f [Hz]') ;
  ylabel('"Amplification" [a.u.]') ;
  if (bPlotLogScale)
    set(gca, 'YScale', 'log') ;
  end
  xlim(FreqRangeShow) ;
  title(sprintf('"Amplification" from: %s%s', AudioFileName, AudioFileExt), ...
        'Interpreter','none') ;

%% Save results:

if (bSaveOutput)

  % Create vaiable to save in .mat file
  f = TransferFuncStruct.f ; % [Hz] - frequency points
  Amplification = TransferFuncStruct.Amplification ; % The trasnfer function
  AmplificationStdErr = TransferFuncStruct.AmplificationStdErr ; % Error in transfer func.
%   % (Renamed parameters for compatibility with previous Amplification files.)
%   FreqsRaster = TransferFuncStruct.f ; % [Hz] - frequency points
%   RRaster = TransferFuncStruct.Amplification ; % The trasnfer function
%   RErrRaster = TransferFuncStruct.AmplificationStdErr ; % Error in transfer func.
  
  if (~isfolder(OutputFolder))
    mkdir(OutputFolder) ;
  end
  
%   save(OutpufFile, 'FreqsRaster', 'RErrRaster', 'RRaster') ;
  save(OutpufFile, 'f', 'Amplification', 'AmplificationStdErr') ;
  
  % Tell user where file was saved.
  fprintf(1, 'Transfer function written to: %s\n', OutpufFile) ;
end

  