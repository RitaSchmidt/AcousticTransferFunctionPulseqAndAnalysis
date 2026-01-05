% -------------------------------------------------------------------------
% GenerateAcousticTrasnferFunctionScan
% -------------------------------------------------------------------------
%
% A script to generate a pulseq .seq file for measuring acoustic (or
% vibrational) trasnfer functions.
%
% Uses a set of triangle shaped gradients to probe the acoustical and/or
% vibrational response of the MRI:
% * The triangle duration defines the BW measured (the response has a null
%   at f = 1/tRamp, where tRamp is half the duration of the triangle).
% * The triangles are played with a constant separation between them,
%   controlled by the user. (Originally a time where they hardly interact
%   with each other, but probably not neccessary.)
% * A "TR" of such triangles is determined automatically from the final
%   desired resolution (in Hz) of the transfer function. (Idealy, 1/TR =
%   resoltuion in Hz.) At the end we will Fourier trasnform (FT) the
%   recorded signal from such "TR"s.
% * For averaging, the "TR" is repeated several times.
% * Because the triangles are repeated within a "TR" the Fourier transform
%   of the recorded signal has the desired frequency resolution, but with
%   many zeros in it. This is because the repeated triangles are the same
%   as a convolution of a single trianlge with a comb of delta function.
%   The FT of this is the same as the FT of a single triangle times(!) a
%   different comb of deltas.
%   To get "in between" data we could add a linear phase to the comb we
%   convolve with, this would shift the position of the final comb (after
%   FT). However, we cannot add a linear phasem because the gradients are
%   real quantities. Instead we apply a sinusoidal envelope. This is
%   effectively the same the as the combination of two linear phases (one
%   positive and one negative).
%
% See Inputs section below for the different parameters and their meaning.
%
% An analysis of the recorded audio/vibrations based on FT of each "TR"
% will extract the transfer function.
%
% For ease of analyis, relevant parameters are written to the [DEFINITIONS]
% section of the .seq file. The requested values are written as comments to
% the file, while the actual (whether different or not) are written without
% being commented.
%
% A .seq file will be saved only if the variable SeqName is not empty.


%% Inputs

% Mark start of preplimiary setups (and initialization and ...)
tStart = tic ;
fprintf(1, 'Preliminary setup ... ') ;


% Desired transfer function resolution (in Hz)
Setup.TransferFuncResHz = 1 ; % [Hz]
% Desired transfer function BW (in Hz). 
% This is only the positive half (since the negative part is symmetric).
% Defines ramp time of triangular gradient.
Setup.TransferFuncPosBWHz = 5000 ; % [Hz]
% Time between start of two gradients (which excite the acoustics).
% Typically, approximate time for sound to decay (neglecting the gradient
% duration itself)
Setup.TimeBetweenGrads = 50e-3 ; % [s]

% Averaging (to improve SNR)
Setup.NumAverages = 10 ; 

% Physical axis
Setup.PhysicalAxis = 'x' ; % 'x', 'y', or 'z'

% NOTE: Up to a linear phase the spectrum of a single triangular gradient
%       of amplitude G and ramp time tRamp is:
%       G*tRamp*(sinc(2*pi*f*tRamp/2)).^2,
%       where f is the frequency.
%       Thus the first zero is at f = 1/tRamp and the BW covered by the
%       center lobe is BW = 2/tRamp (including the negative frequencies
%       which are symmetric, up to a phase).

% Name of sequence, also used as the name to save the .seq file.
% If empty or missing, no .seq file will be saved.
SeqName = 'ScanAcousticTransferFunc' ;

% Should Pulseq plot the resulting sequence?
bSupressPulseqPlots = false ;

%% Initialize Actual Params

% Some parameters used may be different than those requested in the Setup
% structure. These actual parameters will be stored in the structure
% Actual.

% Initialize Actual to be the same as Setup
Actual = Setup ;

%% Read system parameters

% -------------
% Set system HW
% -------------

% Select HW system
[SystemNormal, ...
 SystemsStructArray, ...
 SystemNormalIdx] = Systems.SiemensTerraXR() ; % Terra (7T) defintions


% Set gradient limits for different elements of the sequence

% Default system to use
SystemDefault = SystemNormal ;

GradientRaster = SystemNormal.gradRasterTime ; % [s]


%% Define scan params

% Define single triangular gradient.
% ----------------------------------

tRamp = 1/Actual.TransferFuncPosBWHz ; % [s] Desired!
% round up to raster time
tRamp = ceil(tRamp/GradientRaster) * GradientRaster ;
tPlateau = 0 ; % [s] - For now we set it to zero.
% gradient amplitude (uses maximal slew rate to not exceed maximum allowed
% amplitude).
GradAmp = min(SystemDefault.maxGrad, tRamp*SystemDefault.maxSlew); % [Hz/m]

% Update TransferFuncPosBWHz
Actual.TransferFuncPosBWHz = 1/tRamp ; % [Hz]

% Number of identical repetions - to achieve desired freq. resolution
% -------------------------------------------------------------------

% Number of repetitions of the gradients to run (and FFT) together, to
% reach desired resolution of TransferFuncResHz (or slightly better). Time
% from start of one gradient, to the next, is given by TimeBetweenGrads.
NumGradsPerTR = ceil((1/Actual.TransferFuncResHz)/Actual.TimeBetweenGrads) ;

% Sine-amplitude-modulations
% ---------------------------------

% Number of times to repeat the NumGradsPerTR gradients, but each time with a
% different amplitude modulation.
% Because of TimeBetweenGrads the acoustics (or FFT of the gradient) will
% have many zeros, except every 1/TimeBetweenGrads. In order to shift the
% position of the non-zero samples we could add a linear phase to the
% gradient, but a gradient cannot be complex. This is why we cosine
% amplitude modulate it, so it will be the sum of two complex linear phases
% (each one in a different direction).
NumSineModulations = NumGradsPerTR ;

% The amplitude modulations to use
% NOTE: Once we take the real part there is a redundancy
AmplitudeModulations = real(dftmtx(NumSineModulations)) ;

% "Grad. fill"
% ------------

% Ensure we are on the gradient raster.
Actual.TimeBetweenGrads = ceil(Actual.TimeBetweenGrads / GradientRaster) * ...
                          GradientRaster ;

% Time from end of one gradient to the start of the next to achieve
% TimeBetweenGrads (time between start of two gradients)
GradFill = Actual.TimeBetweenGrads - (2*tRamp +tPlateau) ; % [s]

% Define "TR" 
% -----------

% The TR is not really used, just for reporting.
% When averaging, this is the period of time repeated/averaged.
TR = NumGradsPerTR * Actual.TimeBetweenGrads ; % [s]

% Update TransferFuncResHz
Actual.TransferFuncResHz = 1/TR ; % [Hz]

%% Sanity checks

% NOTE: We assume here that all gradient modes have the same resonance
%       frequecies (and bandwidths).

% Make sure resonance frequencies are defined:
if (~isfield(SystemNormal, 'resonanceFrequencies') || ...
    ~isfield(SystemNormal, 'resonanceBandwidths'))
  error('Resonance frequencies are not defined, you may damage the system') ;
end

% Check we have the same number of resonance frequencies and
% resonance bandwidths.
if (numel(SystemNormal.resonanceFrequencies) ~= ...
                                 numel(SystemNormal.resonanceBandwidths))
  erorr(['Number of resonance frequencies and resonance bandwidth do ' ...
         'not match']) ;
end

% Check we are not applying a forbidden frequency.
% (We just check that the time between triangle gradients does 
NumResonanceFreqs = numel(SystemNormal.resonanceFrequencies) ;
NominalFrequency = 1/Actual.TimeBetweenGrads ; % [Hz]
for ResonanceCounter = 1:NumResonanceFreqs
  ResonanceFreq = abs(SystemNormal.resonanceFrequencies(ResonanceCounter)) ;
  ResonanceBW = abs(SystemNormal.resonanceBandwidths(ResonanceCounter)) ;
  FreqLowBound = ResonanceFreq - ResonanceBW/2 ;
  FreqHighBound = ResonanceFreq - ResonanceBW/2 ;

  if ( (NominalFrequency >= FreqLowBound) && ...
       (NominalFrequency <= FreqHighBound) )
    error(['One over time between gradients (1/%g = %g Hz) falls with a ' ...
           'forbidden frequency range (%g - %g Hz). This may be dangerous ' ...
           'to the gradients'], ...
           Actual.TimeBetweenGrads, NominalFrequency, ...
           FreqLowBound, FreqHighBound) ;
  end
end

%% Initialize Sequence

% Moved here in case we need it early.
% (For example set IDs to some events that are used many times, without
%  change. This speeds up the adding of blocks to the sequence.)
Seq = mr.Sequence(SystemDefault) ;


% Mark end of event defintions
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

%% Event definitions

% Mark start of event defintions
tStart = tic ;
fprintf(1, 'Event definitions ... ') ;

% -------------------------------------------------------------------------
% Basic triangle gradient (before sine-modulation)
% -------------------------------------------------------------------------

eTriangleGrad = mr.makeTrapezoid(Actual.PhysicalAxis, SystemDefault, ...
                                 'amplitude', GradAmp, ...
                                 'riseTime', tRamp, ...
                                 'flatTime', 0, ...
                                 'fallTime', tRamp) ;

% -------------------------------------------------------------------------
% Delay between gradients - "Grad. fill"
% -------------------------------------------------------------------------

eGradFill = mr.makeDelay(GradFill) ; 

% -------------------------------------------------------------------------
% Dummy RF
% -------------------------------------------------------------------------

% IDEA viewer does not seem to visualize anything if there is no RF
% defined, so we define a dummy RF here. Pulseq does not like an RF with
% zero amplitude, at least not a rectangular pulse, so give it a small flip
% angle.

% FlipAngleRad = 1 ; % [radians]
% RFDuration = eTriangleGrad.riseTime + eTriangleGrad.flatTime + ...
%              eTriangleGrad.fallTime + ...
%              -(SystemDefault.rfRingdownTime); % [s]
% [eDummyRF] = mr.makeBlockPulse(FlipAngleRad, SystemDefault, ...
%                                'duration', RFDuration, ...
%                                'use', 'excitation', ...
%                                'delay', SystemDefault.rfDeadTime) ;
FlipAngleRad = 1*pi/180 ; % [radians]
RFDuration = eTriangleGrad.riseTime + eTriangleGrad.flatTime + ...
             eTriangleGrad.fallTime ; % [s]
[eDummyRF] = mr.makeBlockPulse(FlipAngleRad, SystemDefault, ...
                               'duration', RFDuration, ...
                               'use', 'excitation') ;



% Mark end of event defintions
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;


%% Create sequence cell/structre

% Mark start of building sequence blocks
tStart = tic ;
fprintf(1, 'Building sequence ... \n') ;

% Add a dummy RF just so the viewer in IDEA will work.
Seq.addBlock(eDummyRF) ;
% Also add a delay, just so there will be some time before the gradients
% start. (Probably not needed.)
Seq.addBlock(eGradFill) ;

% Loop over different sine-modulations
for AmpModCounter = 1:NumSineModulations
  % Repeat a given sine modulation NumAverages times
  for AvgCounter = 1:Actual.NumAverages
    % Go over NumReps of gradients, but modulate the amplitude according to
    % current modulation.
    for Counter1 = 1:NumGradsPerTR
      % update gradient amplitude
      GAmp = GradAmp * AmplitudeModulations(Counter1, AmpModCounter) ;
      eTriangleGrad.amplitude = GAmp ;
      % Add gradient to sequence
      Seq.addBlock(eTriangleGrad) ;
      % Add delay until next gradient
      Seq.addBlock(eGradFill) ;
    end
  end
end


% Mark end of filling sequence blocks
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;


%% Update [DEFINITIONS] of sequence including dummy comments (the setup)

% Mark start updating defintions
tStart = tic ;
fprintf(1, 'Updating [DEFINITIONS] ... ') ;


% Set definitions to be used by interpreter
clear Defs4Interpreter ; % clear from previous runs.
Defs4Interpreter.Name = SeqName ;
% Defs4Interpreter.FOV = Actual.FOV ; % not relevant here
Defs4Interpreter.kSpaceCenterLine = 0 ; % count from 0
Defs4Interpreter.kSpaceCenterPartition = 0 ; % count from 0
% NOTE: We skip TE and TR becuase currently the UI on the Siemens scanner
%       checks whether the TE and TR are within range, but they are usually
%       off (the range is a single value). Thus for now, the user (and the
%       DICOM images) will see/contain the WRONG TE and TR (in general)
% Defs4Interpreter.TE = Actual.TE ; % not relevant here
% Defs4Interpreter.TR = Actual.TR ; % not relevant here
Defs4Interpreter.PhysicalAxis = Actual.PhysicalAxis ;
Defs4Interpreter.tGradRampUp = tRamp ; % [s]
Defs4Interpreter.tGradPlateau = tPlateau ; % [s]
Defs4Interpreter.tGradRampDown = tRamp ; % [s]
Defs4Interpreter.TimeBetweenGrads = Actual.TimeBetweenGrads ;
Defs4Interpreter.NumGradsPerTR = NumGradsPerTR ;
Defs4Interpreter.NumSineModulations = NumSineModulations ;
Defs4Interpreter.RepetitionTime = TR ; % [s] Pass TR without calling it 'TR'
Defs4Interpreter.NumAverages = Actual.NumAverages ; 
Defs4Interpreter.TransferFuncPosBWHz = Actual.TransferFuncPosBWHz ; % [Hz]
Defs4Interpreter.TransferFuncResHz = Actual.TransferFuncResHz ; % [Hz]

% Add Setup and Actual structs as comments in Seq (and then in .seq file)
% by exploiting the defintions mechanism. In the comments only Actual
% valies which are different from the Setup, will be shown.
% as well as add definition (Defs4Interpreter above) for use by
% interpreter, if it can read them.
Seq = Services.UpdateSeqDefs(Seq, Setup, Actual, Defs4Interpreter) ;


% Mark end of updating defintions
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

%% Save sequence as a .seq file, if a file name/path is defined

if (exist('SeqName', 'var') && ~isempty(SeqName))

  % Mark start writing sequence to file
  tStart = tic ;
  fprintf(1, 'Writing .seq file ... ') ;

  if (exist('SeqFolder', 'var') )
    % ensure SeqFolder exists
    if (~isfolder(SeqFolder))
      mkdir(SeqFolder) ;
    end
    Seq.write(fullfile(SeqFolder, [SeqName, '.seq'])) ;
  else
    Seq.write([SeqName, '.seq']) ;
  end

  % Mark end of writing sequence to file
  tEnd = toc(tStart) ;
  fprintf(1, 'Done. (%g s)\n', tEnd) ;
  
end

%% check whether the timing of the sequence is correct

% Mark start checking sequence timing
tStart = tic ;
fprintf(1, 'Checking sequence timing ... ');

[ok, error_report]=Seq.checkTiming;

% Mark end of writing sequence to file
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

if (ok)
  fprintf(1, 'Timing check passed successfully.\n');
else
  fprintf(1, '\n');
  fprintf(1, 'Timing check failed! Error listing follows:\n');
  fprintf(1, '\n');
  fprintf(1, [error_report{:}]);
  % When  the list is long the heading (that these are errors) is lost so
  % we repeat it at the end
  fprintf(1, '\n');
  fprintf(1, 'Timing check failed! Error listing above this line.\n');
end


%% plot sequence?

if (~bSupressPulseqPlots)

  % Mark start plotting
  tStart = tic ;
  fprintf(1, 'Plotting sequence ... ');
  
  % Plot waveforms/diagram of sequence
  % ----------------------------------
  
% Plot all
Seq.plot() ;   
  
  % Plot k-space trajectory
  % -----------------------
  [ktraj_adc, t_adc, ktraj, t_ktraj, ...
   t_excitation, t_refocusing] = Seq.calculateKspacePP() ;
  figure; 
    plot(t_ktraj.', ktraj.') ;
    hold on ;
  
    xlabel('t [s]') ;
    ylabel('k [1/m]') ;
    title('k-space components as functions of time') ;
  
    % calculateKspacePP() should return physical x, y, z and not logical RO
    % PE and 3D.
    legend('k_x', 'k_y', 'k_z') ;
  
  
%   % Plot k-space trajectory (3D)
%   % ---------------------------
%   figure; 
%     plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:),'b') ;
%     hold on ; 
%     plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.') ; 
%   
%     % calculateKspacePP() should return physical x, y, z and not logical RO
%     % PE and 3D.
%     xlabel('k_x [1/m]') ;
%     ylabel('k_y [1/m]') ;
%     zlabel('k_z [1/m]') ;
%     grid on ;
%     title('3D k-space') ;
%     legend('trajectory', 'ADC')
  
  % Mark end of plotting
  tEnd = toc(tStart) ;
  fprintf(1, 'Done. (%g s)\n', tEnd) ;

end

%% THE END
return ; 


%% TESTING - gradients and resulting spectrum

tRamp = 0.2e-3 ; % [s]
tPlateau = 0 ; % [s]
tESP = 2*tRamp + tPlateau ; % [s]
ETL = 1 ;
DeltaTE = 50e-3 ; % [s]
NumTEs = 1 ;
TSlice = 0 ; % [s]
NumSlices = 1 ;
TR = DeltaTE * NumTEs ; % [s]
t0_slice = 0 ; % [s]
t0_echo = 0 ; % [s]
GradRaster = 10e-6 ; % [s]

[G0] = BuildGrad(tRamp, tPlateau, tESP, ETL, DeltaTE, NumTEs, ...
                  TSlice, NumSlices, TR, t0_slice, t0_echo, GradRaster) ;


% Add "phase"
HzRes = 1 ; % [Hz]
T = 1/HzRes ; % [s]
Reps = round(T/TR) ;

G = zeros(numel(G0)*Reps, Reps) ;


% for Counter = 1:FFTCenterIndex(Reps)
for Counter = 1:Reps

  % Amps = ones(1, Reps) ;
  % Amps(2:2:end) = -1 ;
  Amps = cos(2*pi/Reps*(Counter-1)*(0:(Reps-1))) ;
  
  G(:, Counter) = reshape(G0(:) * Amps, [], 1) ;
end

NumSamples = size(G,1) ;
t = (0:(NumSamples - 1)) * GradRaster ; % [s]



figure ; 
  plot(t, G, 'LineWidth', 1) ;
  xlabel ('t [s]')
  ylabel ('G [a.u.]')  

S = fftshift(fft(G, [], 1 ), 1) ;
f = ((1:NumSamples) - FFTCenterIndex(NumSamples)) / (NumSamples*GradRaster) ;

figure ; 
  plot(f, abs(S)) ;
  xlabel ('f [Hz]')
  ylabel ('|S| [a.u.]')  
  ylim([0 1.1*max(abs(S(:)))])


figure ; 
  plot(f, sum(abs(S), 2)) ;
  xlabel ('f [Hz]')
  ylabel ('|S| [a.u.]')  
  ylim([0 1.1*max(abs(S(:)))])
