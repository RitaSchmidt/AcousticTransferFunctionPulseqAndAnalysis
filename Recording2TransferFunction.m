function [TransferFuncStruct] = Recording2TransferFunction(AudioFilePath, ...
                                                           PulseqInputs, ...
                                                           TR2AcousticSteadyState, ...
                                                           TransferFuncResHz, ...
                                                           TransferFuncPosBWHz, ...
                                                           AudioFormat, ...
                                                           CSVNumHeaderLines, ...
                                                           bDecimalPointCSV)
% [TransferFuncStruct] = Recording2TransferFunction(AudioFilePath, ...
%                                                   PulseqInputs, ...
%                                                   TR2AcousticSteadyState, ...
%                                                   TransferFuncResHz, ...
%                                                   TransferFuncPosBWHz, ...
%                                                   AudioFormat, ...
%                                                   CSVNumHeaderLines)
%
% Determines the transfer function (amplification) of gradients from an
% audio recording of a specialized sequence. 

% The sequence is supposed to be composed of identical (up to amplitude)
% gradient lobes, every TimeBetweenGrads seconds. Every NumGradsPerTR such
% lobes are considered a TR, whos audio recording is FFTed to get an
% acoustic spectrum. Because of the dead time between the gradients, the
% resulting spectrum also has (practically) dead frequencies. To excite the
% "dead frequencies" the TR is repeated with a sine modulation of the
% gradient amplitudes within a TR. (This is a sum of two comples linear
% phase modulations, one shifting the frequency peaks up and one shifting
% the down.) NumSineModulations different sine-modulations are used to
% reach the desired TransferFuncResHz resolution of the trasnfer function.
% It is assumed, for the analysis, that NumGradsPerTR = NumSineModulations.
% For averaging, each TR is repeated NumAverages total times, before moving
% to the next sine modulation.
%
% In the analysis the first TR2AcousticSteadyState "averages" of each TR 
% are discarded, in the assumption that so many TRs are needed for the
% acoustics to reach a steady state (the previous sine-modulation no longer
% affects the current one).
%
% Inputs:
% -------
% - AudioFilePath - Path to audio recording file. Unless AudioFormat is
%   passed, will try to guess file type. (See also  AudioFormat.)
% - PulseqInputs - Inputs on how to analyze the audio recordings. Either a
%   structure with the following fields, or a file (pulse .seq file) with
%   the following keys under the [Definitions] section:
%   > TransferFuncResHz - The frequency resolution of the output transfer
%     function.
%     Can be overridden by an explicit value TransferFuncResHz passed to
%     the function. In such a case the value is not read from the
%     structure/file, so cam be missing as well.
%     Can be missing, if value is given as inpout to function.
%   > TransferFuncPosBWHz - The positive half of the transfer function's
%     BW. (The output trasnfer function will cover frequencies from
%     -TransferFuncPosBWHz to +TransferFuncPosBWHz, up to rounding errors
%     due to the resolution used.) 
%     Can be overridden by an explicit value TransferFuncPosBWHz passed to
%     the function. In such a case the value is not read from the
%     structure/file, so cam be missing as well.
%     Can be missing, if value is given as inpout to function.
%   > PhysicalAxis - A single character 'x', 'y', or 'z', marking which
%     physical gradient was used to generate the audio recording.
%     (The input is not used, it is just passed on to the output
%     structure.)
%   > TimeBetweenGrads - The time in seconds from the begining of one
%     gradient lobe to the following one in the scan.
%   > GradientRaster - The temporal resolution, in seconds, of the
%     gradients on the system. 
%   > NumSineModulations - Number of different sine-modulations used to
%     modulate the gradients within a TR.
%   > NumAverages - Number of time each TR, with the same sine-modulation,
%     is repeated.
%   > NumGradsPerTR - Number of gradient lobes (and the following
%     delay/fill) within a single TR.
%   > RepetitionTime - The TR in seconds. The name RepetitionTime is used
%     instead of TR, to avoid possible UI problems in Pulseq if 'TR' is
%     used.
%   > tGradRampUp - Ramp up duration, in secods, of each trapezoid gradient
%     lobe. 
%   > tGradPlateau - Duration, in seconds, of constant gradient within
%     each trapezoid gradient lobe. (If zero, the shape is a triangle).
%   > tGradRampDown - Ramp up duration, in secods, of each trapezoid 
%     gradient lobe. (Typically, the same as tGradRampUp.)
% - TR2AcousticSteadyState - In the analysis of the audio recording, how
%   many averages of each sine-modulation to discard, so that (assumingly)
%   the acoustics is in steady-state, i.e., not affected by previuos
%   sine-modulation.
%   If after discarding, only a single average remains, the output will
%   have an error estimation of transfer function.
%   [Default: 1] 
% - TransferFuncResHz - An optional value to override the frequency
%   resolution of the output transfer function. If given, a value is not
%   required to be included in PulseqInputs.
%   [Default: [] = no override]
% - TransferFuncPosBWHz - An optional value to override the frequency
%   BW coverage of the output transfer function. If given, a value is not
%   required to be included in PulseqInputs. The output transfer function
%   will cover the frequencie -TransferFuncPosBWHz to +TransferFuncPosBWHz
%   up to the accuracy of frequency resolution used.
%   [Default: [] = no override]
% - AudioFormat - A string matching the format of the audio file. Can be
%   one of the following (case insensitive):
%   > 'AudioFile' - An audio file format recognized by audioread().
%   > 'TimeTable' - A .mat file expected to include only a single timetable
%     variable. The timetable is expected to include only two columns, one
%     of which is called 'time' (case insensitive). The timetable is also
%     expected to have a property 'SampleRate' (case sensitive!) in Hz.
%   > 'CSV' - A CSV formatted file, with a header and exactly two columns,
%     the first of which is considered to be time in seconds.
%     The header text, CSVNumHeaderLines lines long, is ignored.
% - CSVNumHeaderLines - Number of lines to skip at the start of a CSV
%   format audio file. If not given, or empty, a default is used (16).
%   [Default: [] = used default of 16]
% - bDecimalPointCSV - In case of a CSV audio file, determines decimal
%   separator used. If true, assumes a decimal point (e.g., 1+1/10 = 1.1).
%   If false assumes a decimal comma (1.1 --> 1,1). If empty(!) will try to
%   determine format on its own.
%   [Default: [] = try to guess decimal separator]

%
% Outputs:
% --------
% - TransferFuncStruct - A structre containing the resulting trasnfer
%   function. It includes the following fields
%   > f - The frequencies, in Hz, at which the transfer function is
%     defined.
%   > Amplification - The actual transfer function. A relative
%     amplification due to the mechanical resonances.
%   > AmplificationStdErr - An estimated standard error for the transfer
%     function amplification. If there are not enough averages (after
%     discarding TR2AcousticSteadyState averages, this will be empty.
%   > PhysicalAxis - A string stating to which physical axis the transfer
%     function applies. Should be 'x', 'y', or 'z', but there is no testing
%     of this.
  
  % To do:
  % ------
  % * Fix linear phase of each TR (per sine-modulation) based on non-zero
  %   points in the model (or practically non-zero).
  % * Could we use the magnitude of the spectrum std as an additional
  %   measure of the amplification. (The FFT noise should be uniform, so
  %   any difference in std is due to amplification.)
  % * Properly normalize!
  % * Add help to sub-functions.
  

%% Setup

switch nargin
  case 2
    TR2AcousticSteadyState = [] ; % Do not override value from PulseqInputs
    TransferFuncResHz = [] ;      % Do not override value from PulseqInputs
    TransferFuncPosBWHz = [] ;    % Do not override value from PulseqInputs
    AudioFormat = [] ;            % Guess by extension
    CSVNumHeaderLines = [] ;      % Use default set in ReadAudio(), if needed.
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 3
    TransferFuncResHz = [] ;      % Do not override value from PulseqInputs
    TransferFuncPosBWHz = [] ;    % Do not override value from PulseqInputs
    AudioFormat = [] ;            % Guess by extension
    CSVNumHeaderLines = [] ;      % Use default set in ReadAudio(), if needed.
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 4
    TransferFuncPosBWHz = [] ;    % Do not override value from PulseqInputs
    AudioFormat = [] ;            % Guess by extension
    CSVNumHeaderLines = [] ;      % Use default set in ReadAudio(), if needed.
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 5
    AudioFormat = [] ;            % Guess by extension
    CSVNumHeaderLines = [] ;      % Use default set in ReadAudio(), if needed.
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 6
    CSVNumHeaderLines = [] ;      % Use default set in ReadAudio(), if needed.
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 7
    bDecimalPointCSV = [] ;       % Try to automatically detect (if empty).
  case 8
    % Do nothing.
  otherwise
    error('Too many or too few input arguments.') ;
end

% BaseFolder = '/home/amir/Academic/Data/ForbiddenFrequenciesAndTimingWithinTR' ;
% RecordingFolder = fullfile(BaseFolder, '2025.06.05-Pulseq-Amplification-Test', ...
%                            'Recordings') ;
% 
% % RecordingFile = fullfile(RecordingFolder, 'MeasUID-1029.wav') ; % 2 averages
% RecordingFile = fullfile(RecordingFolder, 'MeasUID-1030.wav') ; % 10 averages
% 
% TransferFuncPosBWHz = 5000 ; % [Hz] - positive half
% TransferFuncResHz = 1 ; % [Hz]
% PhysicalAxis = 'x' ; % 'x', 'y', or 'z'
% TimeBetweenGrads = 0.05 ; % [s]
% NumSineModulations = 20 ;
% % NumAverages = 2 ;
% NumAverages = 10 ;
% NumGradsPerTR = 20 ;
% TR = 1 ; % [s]
% % How many TRs (we assume) it takes to reach an acoustic steady-state. We
% % will remove this number of the initial averages from the signal used.
% TR2AcousticStadState = 1 ;
% 
% tGradRampUp = 0.2e-3 ; % [s]
% tGradPlateau = 0 ; % [s]
% tGradRampDown = 0.2e-3 ; % [s]

% GradientRaster = 10e-6 ; % [s]

  %% Read defintions from PulseqInputs

  if isstruct(PulseqInputs)

    if isempty(TransferFuncPosBWHz)
      TransferFuncPosBWHz = PulseqInputs.TransferFuncPosBWHz ; % [Hz]
    end
    if isempty(TransferFuncResHz)
      TransferFuncResHz   = PulseqInputs.TransferFuncResHz ; % [Hz]
    end
    PhysicalAxis        = PulseqInputs.PhysicalAxis ; % 'x', 'y', or 'z'
    TimeBetweenGrads    = PulseqInputs.TimeBetweenGrads ; % [s]
    GradientRaster      = PulseqInputs.GradientRasterTime ; % [s]
  
    NumSineModulations  = PulseqInputs.NumSineModulations ;
    NumAverages         = PulseqInputs.NumAverages ;
    NumGradsPerTR       = PulseqInputs.NumGradsPerTR ;
    TR                  = PulseqInputs.RepetitionTime ; % [s]
    tGradRampUp         = PulseqInputs.tGradRampUp ; % [s]
    tGradPlateau        = PulseqInputs.tGradPlateau ; % [s]
    tGradRampDown       = PulseqInputs.tGradRampDown ; % [s]

  else % read data from a pulseq file
    % Read defintions section from .seq file
    Defs = DefsFromSeq(PulseqInputs) ;
  
    if isempty(TransferFuncPosBWHz)
      TransferFuncPosBWHz = Defs.TransferFuncPosBWHz ; % [Hz]
    end
    if isempty(TransferFuncResHz)
      TransferFuncResHz = Defs.TransferFuncResHz ; % [Hz]
    end
    PhysicalAxis = Defs.PhysicalAxis ; % 'x', 'y', or 'z'
    TimeBetweenGrads = Defs.TimeBetweenGrads ; % [s]
  
    GradientRaster = Defs.GradientRasterTime ; % [s]

    NumSineModulations = Defs.NumSineModulations ;
    NumGradsPerTR = Defs.NumGradsPerTR ;
    TR = Defs.RepetitionTime ; % [s]
    NumAverages = Defs.NumAverages ; 
    tGradRampUp = Defs.tGradRampUp ; % [s]
    tGradPlateau = Defs.tGradPlateau ; % [s]
    tGradRampDown = Defs.tGradRampDown ; % [s]
    

%     warning('Some inputs are hardcoded!!!!!')
%   
%     NumSineModulations = 20 ;
%     NumAverages = 10 ;
%     NumGradsPerTR = 20 ;
%     TR = 1 ; % [s]
%     tGradRampUp = 0.2e-3 ; % [s]
%     tGradPlateau = 0 ; % [s]
%     tGradRampDown = 0.2e-3 ; % [s]
  end




  % How many TRs (we assume) it takes to reach an acoustic steady-state. We
  % will remove this number of the initial averages from the signal used.
  if isempty(TR2AcousticSteadyState)
    TR2AcousticSteadyState = 1 ;
  end
  

  %% Sanity checks

  % Are NumSineModulations and NumGradsPerTR not the same?
  if (NumSineModulations ~= NumGradsPerTR)
    error('NumSineModulations (%d) is not the same as NumGradsPerTR (%d).', ...
          NumSineModulations, NumGradsPerTR) ;
  end


  %% read audio file and split it into TRs
  
  % if isempty(CSVNumHeaderLines)
  %   [Audio, Fs] = ReadAudio(AudioFilePath, AudioFormat) ;
  % else
  %   [Audio, Fs] = ReadAudio(AudioFilePath, AudioFormat, ...
  %                           CSVNumHeaderLines, bDecimalPointCSV) ;
  % end
  [Audio, Fs] = ReadAudio(AudioFilePath, AudioFormat, ...
                          CSVNumHeaderLines, bDecimalPointCSV) ;
  
  dt = 1/Fs ; % [s]

  % Remove quiet periods before/after the scan in the recording and split
  % it into TRs
  [AudioSplit] = TrimAndSplitAudio2TRs(Audio, dt, TimeBetweenGrads, ...
                                                NumSineModulations, ...
                                                NumAverages, TR) ;  
  



  % Number of audio samples per TR. In general the number might not be an
  % integer, which will result in a slight shift between audio and actual
  % scan. Thus we find the number of audio samples per TR twice. A floating
  % number (non-integer) and a rounded value
  SamplesPerTR_NonInt = TR/dt ;
  SamplesPerTR = round(SamplesPerTR_NonInt) ;
  
  % Remove "dummy" intial averages of each sine-modulation, so the audio
  % used will be in acoustic steady-state (not affected from the previous
  % sine-modulation).
  % For now the number of initial averages is given, we could(?) do
  % something smarter, based on the recordings themselves/
  % Number of averages actually used
  if islogical(TR2AcousticSteadyState)
    bAveragesUse = false(NumAverages, 1) ;
    NumAveragesModify = min(NumAverages, numel(TR2AcousticSteadyState)) ;
    bAveragesUse(1:NumAveragesModify) = ...
                              TR2AcousticSteadyState(1:NumAveragesModify) ;
    AudioSplit = AudioSplit(:,bAveragesUse, :) ;
    NumNetAverages = nnz(bAveragesUse) ;
    
  else
    AudioSplit = AudioSplit(:,(TR2AcousticSteadyState+1):end, :) ;
    NumNetAverages = NumAverages - TR2AcousticSteadyState ;
  end

  %% Generate measured spectra


  SpectrumPerTR = fftshift(fft(AudioSplit, [], 1), 1) ;
  CenterSpecrtrumIdx = FFTCenterIndex(SamplesPerTR) ; % counting from 1
  % NOTE: dt*SamplesPerTR might not be exactly a TR, because (SamplesPerTR 
  %       is forced to be an integer).
  f = ((1:SamplesPerTR) - CenterSpecrtrumIdx) / (dt*SamplesPerTR) ;


  %% Initial spectrum linear phase correction - due to SamplesPerTR_NonInt
  
  % We will only (try to) fix the shift between the aveages of each
  % sine-modulation, because we will use the averages to find a mean and
  % standard deviation per sine-modulation. The different sine-modulations
  % will be combined in their absolute values, so we don't care(?) about
  % their relative phases.
  
  IdxShiftPerTR = SamplesPerTR - SamplesPerTR_NonInt ;
  
  if (abs(IdxShiftPerTR) > 2*eps(SamplesPerTR))
    % Define phase correction per average within each sine modulation
    PhaseCorr = (2*pi*IdxShiftPerTR/SamplesPerTR) * ...
                (0:(NumNetAverages-1)) .* ...
                ((1:SamplesPerTR) - CenterSpecrtrumIdx).' ;
    % apply the correction
    SpectrumPerTR = SpectrumPerTR .* exp(-1i*PhaseCorr) ;
  end

  %% Limit measured spectrum and frequency to relevant frequencies only

  % The gradient blips in the scan excite a sinc like spectrum
  % (disregarding the mechanical resonances). At the nulls, the
  % amplification found here will not be accurate, so we take only the
  % center lobe of the spectrum's sinc distribution.
  
  % Limit to relevant frequencies only
  bInSupportedFreq = (f >= -abs(TransferFuncPosBWHz) & ...
                      f <= abs(TransferFuncPosBWHz)) ;
  f = f(bInSupportedFreq) ;
  SpectrumPerTR = SpectrumPerTR(bInSupportedFreq, :,:) ;

  %% interpolate to desired frequency resolution

  % We could interpolate the final transfer function, but I do not know how
  % to interpolate the errors in the case. Therefore we interpolate the
  % "raw" spectra before any averaging (and calculating std). This way the
  % errors will come out automatically

  % interpolate only if TransferFuncResHz > 0
  if (TransferFuncResHz > 0)
    NFreq = floor(TransferFuncPosBWHz / TransferFuncResHz) ;
    fOut = ((-NFreq):NFreq) * TransferFuncResHz ;
    
    % We want to interpolate multiple columns, so we will use
    % griddedInterpolant.
    InterpMethod = 'pchip' ; % <---- Good enough?
    F = griddedInterpolant(f(:), SpectrumPerTR, InterpMethod) ; 

    % Update SpectrumPerTR and f
    SpectrumPerTR = F(fOut(:)) ;
    f = fOut(:) ;
  end

  %% Build model spectrum (no amplification due to mechanical resonances)
  

  ModelSpectra = GetModelSpectra(f, tGradRampUp, tGradPlateau, ...
                                 tGradRampDown, ...
                                 TimeBetweenGrads, NumGradsPerTR, ...
                                 NumSineModulations, ...
                                 GradientRaster, 'pchip') ;



  %% Average and standard deviation of spectra, per sine modulation
  
  % We take the complex(!) average and after averaging use only the absolute
  % value.
  AvgMeasSpectra = reshape(abs(mean(SpectrumPerTR, 2)), ...
                           numel(f), NumSineModulations) ;
  if (NumNetAverages > 1)
    StdErrMeasSpectra = reshape(std(SpectrumPerTR, 0, 2), ...
                                numel(f), NumSineModulations) / ...
                        sqrt(NumNetAverages) ;
  else
    StdErrMeasSpectra = [] ; 
  end
  
  % Define mask (per sine-modulation) of which frequencies to use for a
  % linear phase correction
  MaxAbsModel = max(abs(ModelSpectra(:))) ;
  bValidFreqPerSineModulation = abs(ModelSpectra) > 10*eps(MaxAbsModel) ;

  % Find, unlikely case, of std also being zero (when there are averages)
  % In this case, for he given frequency we will use only the model for
  % weighting.
  if (NumNetAverages > 1)
    bValidStd1D = any(StdErrMeasSpectra> 10*eps(MaxAbsModel), 2) ; 
    bValidStd2D = repmat(bValidStd1D, 1, NumSineModulations) ;
  end


  %% Weighted average of amplification
  
  % To avoid devision by zero, we define the inverse of the standard
  % deviation and the model. Where the model is close enough to zero (not
  % bValidFreqPerSineModulation), the inverses are defined as zero. (As a
  % result they won't contribute to the calculation.)
  if (NumNetAverages > 1)
    InvStdErr = zeros(size(AvgMeasSpectra)) ;
    bValidFreqsAndStd = bValidFreqPerSineModulation & bValidStd2D ;
    InvStdErr(bValidFreqsAndStd) = 1./StdErrMeasSpectra(bValidFreqsAndStd) ;
  end
  InvModel = zeros(size(AvgMeasSpectra)) ;
  InvModel(bValidFreqPerSineModulation) = ...
                               1./ModelSpectra(bValidFreqPerSineModulation) ;
  
  % The amplification is the ratio between the measured signal and the model.
  % We weight this with the inverse(!) of expected standard deviation of the
  % ratio (of the averages per sine-modulation). It will be used only
  % squared, so we square it already now.
  % If, however, there not enough averages (i.e., just one real
  % measurement), the standard deviation is meaningless, so don't use it.
  if (NumNetAverages > 1)
    WeightSqrd = zeros(numel(f), NumSineModulations) ;
    WeightSqrd(bValidStd1D, :) = abs(InvStdErr(bValidStd1D,:) .* ...
                                     ModelSpectra(bValidStd1D,:)).^2 ;
    WeightSqrd(~bValidStd1D, :) = abs(ModelSpectra(~bValidStd1D,:)).^2 ;
    
  else
    WeightSqrd = abs(ModelSpectra).^2 ;
  end
  
  
  Amplification = sum(abs(AvgMeasSpectra .* InvModel) .* WeightSqrd, 2) ./ ...
                  sum(WeightSqrd, 2) ;
  if (NumNetAverages > 1)
    AmplificationStdErr = sqrt(1 ./ sum(WeightSqrd, 2)) ;
    AmplificationStdErr(~bValidStd1D) = NaN ; % no error defined.
  else
    AmplificationStdErr = [] ;
  end


  %% Place results in output structure
  
  TransferFuncStruct.f = f ; % [Hz]
  TransferFuncStruct.Amplification = Amplification ; 
  TransferFuncStruct.AmplificationStdErr = AmplificationStdErr ;
  TransferFuncStruct.PhysicalAxis = PhysicalAxis ;



  %% The END
  return ;

end

%% Helper functions

% Read audio file
function [Audio, Fs] = ReadAudio(AudioFilePath, AudioFormat, ...
                                 CSVNumHeaderLines, bDecimalPointCSV)

  % Default AudioFormat. Empty means, guess by extension.
  AudioFormatDefault = [] ; 

  % Default number of lines at start of CSV file which serve as a header
  % and should be ignored. Can be overridden, if desired.
  CSVNumHeaderLinesDefault = 16 ;
  
  % Default assumption on format on CSV file (if used). If true assumes a
  % decimal point is used in the numbers, if false, assumes a decimal comma
  % (European/non-USA style). If empty(!), will try to determine the format
  % automatically.
  bDecimalPointCSVDefault = [] ; % empty means try to guess.


  % Make sure all inputs exist. Missing inputs are set to empty which is
  % replaced by default values after the switch statement.
  switch nargin
    case 1
      AudioFormat = [] ;
      CSVNumHeaderLines = [] ; 
      bDecimalPointCSV = [] ;
    case 2
      CSVNumHeaderLines = [] ; 
      bDecimalPointCSV = [] ;
    case 3
      bDecimalPointCSV = [] ;
    case 4
      % nothing to do.
    otherwise
      error('Too many or too few input arguments.')
  end

  % Replace empty value with default ones
  if isempty(AudioFormat)
    AudioFormat = AudioFormatDefault ;
  end
  if isempty(CSVNumHeaderLines)
    CSVNumHeaderLines = CSVNumHeaderLinesDefault ;
  end
  if isempty(bDecimalPointCSV)
    bDecimalPointCSV = bDecimalPointCSVDefault ;
  end

  % =======================================================================
  % Find audio format, if not given
  % =======================================================================
  if (isempty(AudioFormat))
    % Guess format by extension
    [~, ~, AudioExt] = fileparts(AudioFilePath) ;

    % Currently(?) supported formats by audioread().
    AudioReadExt = {'.aifc', '.aiff', '.aif', '.au', '.flac', '.ogg', ...
                    '.opus', '.wav', '.mp3', '.m4a', '.mp4'} ;


    switch (lower(AudioExt))
      case lower(AudioReadExt)
        AudioFormat = 'AudioFile' ;
      case lower('.mat')
        AudioFormat = 'TimeTable' ;
      case lower('.csv')
        AudioFormat = 'CSV' ;
      otherwise
        error(['Could not guess audio file type from extension, pleas pass ' ...
               'audio format explicitly as an input']) ;
    end
  end

  % =======================================================================
  % Actual reading of audio file
  % =======================================================================

  switch lower(AudioFormat)
    case lower('AudioFile')
      [Audio, Fs] = audioread(AudioFilePath) ;
    case lower('TimeTable')
      % warning('Reading audio from Matlab table file') ;
      RecordingTable = load(AudioFilePath) ;
      TableName = fieldnames(RecordingTable) ;
      % We expect a single timetable variable to be in the file
      if ~(numel(TableName) == 1 && ...
           isa(RecordingTable.(TableName{1}), 'timetable' ))
        error('Content of file is not a single timetable: %s', ...
              AudioFilePath) ;
      end
      RecordingTable = RecordingTable.(TableName{1}) ;
      TableColumns = RecordingTable.Properties.DimensionNames ;
      % We expecte two dimensions, one of which is Time.
      if (numel(TableColumns) ~= 2)
        error('Time table is expected to include exactly two columns.')
      end
      bTimeColumn = strcmpi(TableColumns, 'time') ;
      % TimeColStr = TableColumns{bTimeColumn} ;
      AudioColStr = TableColumns{~bTimeColumn} ;
      Audio = double(RecordingTable.(AudioColStr)) ;
      Fs = RecordingTable.Properties.SampleRate ;
    case lower('CSV')
      % % % assumes two columns in the CSV, the first is time in seconds and the
      % % % second is the sound level.
      % % Audio = readmatrix(AudioFilePath, 'NumHeaderLines', CSVNumHeaderLines) ;
      % % % We expect an Nx2 matrix to be read
      % % if ~(ismatrix(Audio) && size(Audio, 2) == 2)
      % %   error('CSV file is expected to include an Nx2 matrix')
      % % end
      % % Fs = 1/diff(Audio(1:2,1)) ; % [Hz]
      % % Audio = Audio(:,2) ; % Keep only second column.
      if (isempty(bDecimalPointCSV))
        % Try to guess which mark is used as a decimal separator. (This is
        % not fool proof.)
        % We open the file and read the first line. If it includes ';', we
        % assume this is the separator between values, instead of a comma,
        % which means the comma is reserved for the decimal separator
        % (instead of a period).
        % open file
        fid = fopen(AudioFilePath, 'rt') ;
        % Read first line
        CSVLine = fgetl(fid) ;
        % Does line contain ';' ==> uses comma is decimal separator
        % Or contains ',' ==> uses period is decimal separator
        if (contains(CSVLine, ';'))
          bUseDecimalPointInCSV = false ;
        elseif (contains(CSVLine, ','))
          bUseDecimalPointInCSV = true ;
        else
          error(['Cannot guess if CSV uses decimal period or decimal ', ...
                 'comma. Please set bUseDecimalPointInCSV appropriately, ' ...
                 'for the file.']) ;
          
        end
  
  
        fclose(fid) ;
  
      else
        bUseDecimalPointInCSV = bDecimalPointCSV ;
      end
  
      % Now read the file according to its format.
      if (bUseDecimalPointInCSV) % use a period as a decimal separator (e.g. 3.14)
        % assumes two columns in the CSV, the first is time in seconds and the
        % second is the sound level.
        Audio = readmatrix(AudioFilePath, ...
                           'NumHeaderLines', CSVNumHeaderLines) ;
        % We expect an Nx2 matrix to be read
        if ~(ismatrix(Audio) && size(Audio, 2) == 2)
          error('CSV file is expected to include an Nx2 matrix')
        end
      else % use a comma as a decimal separator (e.g. 3,14) 
        % Read file
        fid = fopen(AudioFilePath, 'rt') ;
        StrCSV = fread(fid, '*char').' ;
        fclose(fid) ; % close file, no longer needed.
        % Translate CSV from decimal-commma format to decimal-period format.
        % 1. replace all commas (,) with period (.)
        StrCSV = strrep(StrCSV, ',', '.') ;
        % 2. replace all semicolons (;) with commas (,)
        StrCSV = strrep(StrCSV, ';', ',') ;
      
        % Now we can read the content, skipping first CSVNumHeaderLines
        % We assume two value per line. No sanity check!
        Audio = textscan(StrCSV, '%f, %f', 'HeaderLines', CSVNumHeaderLines) ;
        try
          Audio = cell2mat(Audio) ;
        catch
          error('Failed reading audio file as CSV with decimal comma(!).') ;
        end
      end
      Fs = 1/diff(Audio(1:2,1)) ; % [Hz]
      Audio = Audio(:,2) ; % Keep only second column.
  
      
    otherwise
      error('Unsupported audio format ''%s''.',AudioFormat)
  end

end

% Trim quiet time around the scan in the audio + split it into TRs.
function [AudioSplit] = TrimAndSplitAudio2TRs(Audio, dt, TimeBetweenGrads, ...
                                              NumSineModulations, ...
                                              NumAverages, TR)

  %% Find start and end of actual scan in recording

  % The scan includes multiple TRs, in each a gradient is switched on and
  % off every TimeBetweenGrads, so we are going to look for this structure
  % in the recording, or actually it start. We will use convolutions for
  % this. Since a convolution is a heavy operation, we will do it first on
  % a "low resolution" recording.


  % Generate a low resolution "audio" by averaging every
  % TimeBetweenGrads/dt/2 samples of the full audio. (The averaging is over
  % the absolute value of the samples.)
  SamplesPerLowResSample = round(TimeBetweenGrads/dt/2) ;
  NumLowResSamples = floor(numel(Audio)/SamplesPerLowResSample) ;
  NumTrimmedSamples = NumLowResSamples*SamplesPerLowResSample ;
  LowResAbsAudio = mean(reshape(abs(Audio(1:NumTrimmedSamples)), ...
                                SamplesPerLowResSample, []), 1) ;
  % Time of actual scan (without extra recording and the dummy RF block)
  TScan = TR * (NumSineModulations*NumAverages) ;
  % Numer of low-res samples in actual scan
  LowResScanSamples = round(TScan/dt/SamplesPerLowResSample) ;
  
  % convolove a rect of duration matching the scan time with the (absolute
  % of the recroding). This should be maximal when scan duration covers the
  % correct region of the recording. 
  % We start by convolving the Low res case (to get a fast initial guess)
  % and then we can convolve a smaller range of the original recording.
  
  LowResConv = conv(LowResAbsAudio(:), ones(LowResScanSamples, 1),  'valid') ;
  % Find maximum index (where convolution is maximal)
  [~, LowResMaxIdx] = max(LowResConv) ;
  
  % Trim the audio based on the low-res convolution
  StartIdx = max(1, (LowResMaxIdx-2)*SamplesPerLowResSample) ;
  EndIdx = min(numel(Audio), ...
               StartIdx + (LowResScanSamples+2)*SamplesPerLowResSample) ;
  AudioTrimmed = Audio(StartIdx:EndIdx) ;
  
  % Find more precise start and end. (After timming the audio, the retangle
  % we convolve the audio with is now similar in size to the audio itself.
  % Thus, the convolution operation should hopefully not be too heavy.)
  HiResConv = conv(abs(AudioTrimmed(:)), ...
                   ones(round(TScan/dt), 1), 'valid') ;
  
  % Define start of scan in audio recording (when convolution is maximal)
  [~, startScanIdx] = max(HiResConv) ;
  startScanIdx = startScanIdx + StartIdx - 1 ;
  % EndScanIdx = startScanIdx + round(TScan/dt) - 1 ;

  % Number of audio samples per TR. In general the number might not be an
  % integer, which will result in a slight shift between audio and actual
  % scan. Thus we find the number of audio samples per TR twice. A floating
  % number (non-integer) and a rounded value
  SamplesPerTR_NonInt = TR/dt ;
  SamplesPerTR = round(SamplesPerTR_NonInt) ;
  
  % We initially split the audio into NumSineModulations segments. Within
  % each segment we allow the start of the TR to drift because of the
  % difference between SamplesPerTR and SamplesPerTR_NonInt. However, we try
  % to compensate for the drift at the start of every sine modulation.
  AudioSplit = zeros(SamplesPerTR*NumAverages, NumSineModulations) ;
  % Set starting (integer) indices of every sine-modulation which best
  % compensate for the drift due to SamplesPerTR vs. SamplesPerTR_NonInt.
  StartSineModulationIdxs = startScanIdx + ...
                             round((0:(NumSineModulations - 1)) * ...
                                   (SamplesPerTR_NonInt * NumAverages)) ;
  % Now fill AudioSplit
  for SineModCounter = 1:NumSineModulations
    StartIdx = StartSineModulationIdxs(SineModCounter) ;
    EndIdx = StartIdx + SamplesPerTR*NumAverages - 1 ;
    AudioSplit(:, SineModCounter) = Audio(StartIdx:EndIdx) ;
  end
  
  % Now we can split the audio into TRs as well
  AudioSplit = reshape(AudioSplit, SamplesPerTR, NumAverages, ...
                       NumSineModulations) ;
  
end

% Generate model spectra 
% (From gradient waveform - without mechanical resonances)
function [ModelSpectra] = GetModelSpectra(f, tGradRampUp, tGradPlateau, ...
                                          tGradRampDown, ...
                                          TimeBetweenGrads, NumGradsPerTR, ...
                                          NumSineModulations, ...
                                          GradientRaster, InterpMethod)

  %% Build model spectrum (no amplification due to mechanical resonances)
  
  
  % Build single trapezoid.

  % GRampUp = (1:round(tGradRampUp/GradientRaster)) - 0.5 ;
  % GRampDown = (round(tGradRampDown/GradientRaster):-1:1) - 0.5 ;
  % GPlateau = ones(round(tGradPlateau/GradientRaster), 1) ;
  % GradTrap = [GRampUp(:)
  %             GPlateau(:)
  %             GRampDown(:)] ;
  GRampUp = (0:round(tGradRampUp/GradientRaster-2)) ;
  GPlateau = ones(round(tGradPlateau/GradientRaster-1), 1) ;
  GRampDown = (round(tGradRampDown/GradientRaster-1):-1:1) ;
  GradTrap = [GRampUp(:)
              GPlateau(:)
              GRampDown(:)] ;
  
  % Generate trapezoid + delay before next trapezoid
  GSubTR = zeros(round(TimeBetweenGrads/GradientRaster), 1) ;
  GSubTR(1:numel(GradTrap)) = GradTrap ;
  
  % replicate trapezoid + delay to complete a TR and scan + apply sine
  % modulations
  Modulations = dftmtx(NumGradsPerTR) ;
  
  GScan = repmat(GSubTR, 1, NumGradsPerTR, NumSineModulations) .* ...
          reshape(real(Modulations), 1, NumGradsPerTR, NumSineModulations) ;
  GScan = reshape(GScan, [], NumSineModulations) ;
  

  % Generate spectrum
  SpectraTheory = fftshift(fft(GScan, [], 1), 1) ;
  
  LentghSpectrumTheory = size(SpectraTheory, 1) ;
  fTheory = ...
    ((1:LentghSpectrumTheory) - FFTCenterIndex(LentghSpectrumTheory)) / ...
    (GradientRaster * LentghSpectrumTheory) ;
  
  % interpolate to recorded spectrum frequencies
  ModelSpectra = complex(zeros(numel(f), NumSineModulations)) ;
  for SineModCounter = 1:NumSineModulations
    ModelSpectra(:, SineModCounter) = ...
           interp1(fTheory, SpectraTheory(:,SineModCounter), f, ...
                   InterpMethod, 0) ;
  end

end

% Read [Definitions] from .seq file
function [Defintions] = DefsFromSeq(SeqFilePath)

   % Read all lines of file in to a cell array
   [TextLines] = textread(SeqFilePath,'%s','delimiter','') ;
   NumLines = length(TextLines) ;
  
   % Go over and analyze each line
   bSectionDefs = false ;
   
   for LineCounter = 1:NumLines
      
      if ~isempty(regexp(TextLines{LineCounter}, ...
                         '^\s*\[\s*(\S*)\s*\]\s*$', 'once'))
      % line defines a section header
         SectionString = regexp(TextLines{LineCounter}, ...
                                '^\s*\[\s*(\S*)\s*\]\s*$', 'once', ...
                                   'tokens') ;
         if strcmpi(SectionString{1}, 'DEFINITIONS')
           bSectionDefs = true ; 
         else
           if (bSectionDefs)
             % We are starting a new section after we have read the
             % [DEFINITIONS] section, so we completed reading what we
             % looked for.
             return ;
           else
             bSectionDefs = false ; % should not be necessary.
           end
         end
                                                 
      elseif (bSectionDefs && ...
              ~isempty(regexp(TextLines{LineCounter},'^s*(#)', 'once')))
      % line is commented (starts with '#')
         % Do nothing, i.e. skip line because it a comment
         
      elseif (bSectionDefs && ...
              ~isempty(regexp(TextLines{LineCounter}, ...
                             '^s*(\S+)\s*(\S+)', 'once')) )
      % line is a variable-value pair (with a white space between them)
        LineStr = TextLines{LineCounter} ;
        TokenExtents = regexp(LineStr, '^s*(\S+)\s*(\S+)', 'once', ...
                                      'tokenExtents') ;
        StartIdxs = TokenExtents(:,1) ;
        EndIdxs = TokenExtents(:,2) ;
        Key = LineStr(StartIdxs(1):EndIdxs(1)) ;
        % Now are the values set to the key numeric or text (we assume not
        % a combination)
        if (isempty(sscanf(LineStr(StartIdxs(2):EndIdxs(2)), '%g')))
          % It is an array of strings
          Values = textscan(LineStr(StartIdxs(2):end), '%s') ;
          if iscell(Values)
            Values = Values{1} ;
            if (iscell(Values) && numel(iscell(Values)) == 1)
              Values = Values{1} ; % make it a simple string.
            end
          end
        else
          % It is an array of numbers
          Values = sscanf(LineStr(StartIdxs(2):end), '%g') ;
        end

        % Store pair
        Defintions.(Key) = Values ;
      
      else
        if (bSectionDefs) % only a problem if we are in [DEFINITIONS]
         error('Do not know to handle inifile line:\n %s', ...
               TextLines{LineCounter}) ;
        end
      end
   end
return ;  


end





