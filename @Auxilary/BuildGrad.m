function [G] = BuildGrad(tRamp, tPlateau, tESP, ETL, DeltaTE, NumTEs, ...
                         TSlice, NumSlices, TR, t0_slice, t0_echo, GradRaster)
% [G] = BuildGrad(tRamp, tPlateau, tESP, ETL, DeltaTE, NumTEs, ...
%                 TSlice, NumSlices, TR, t0_slice, t0_echo, GradRaster)
%
% Create a gradient waveform of the readout gradient-trains in a single TR
% of EPI. Each echo (TE) of each slice in the TR is represented by a single
% gradient train. Multiple such trains are generated for a given TR.
% All times should be given in the same units. All times are assumed to be
% integer multiples of GradRaster and are rounded to the nearest
% GradRaster.
%
% Inputs
% ------
% - tRamp - Time it takes for the trapezoid gradients to ramp up and ramp
%   down (assumed the same).
% - tPlateau - Time the trapezoid gradient is constant (from the time the
%   ramp up ended to the time the ramp down starts.
% - tESP - Time from the start of one gradient to the start of the next one
%   (of opposite sign). Usually will be the duration of a single trapezoid,
%   but is allowed to be longer.
% - ETL - 'Echo Train Length'. Total number of trapezoids (positive and
%   negative) in the gradient train. ETL*tESP is the duration of the
%   gradient train.
% - DeltaTE - Time between start of one gradient train to the start of the
%   next within a single slice (when there are multiple TEs). Must be
%   given, but used only when NumTEs > 1.
% - NumTEs - Numer of echos, TEs, per slice. Must be one or greater. 
%   See also DeltaTE.
% - TSlice - Time per slice. Time from start of first gradient train to the
%   startof the first gradient train of the next slice. Must be given, but 
%   used only when NumSlices > 1.
% - NumSlices - Number of slices per TR. Must be 1 or greater. (Assumes no
%   simulatenous-multi-slicing (SMS). Instead counts all parallel excited
%   slices as a single slice.)
% - TR - Repetition Time - Time allowed for acquiring all slices. This will
%   be the total time of the gradient waveform generated. (For multiple
%   TRs, just replicate this waveform.)
% - t0_slice - Starting time of the first slice within a TR. )See also
%   t0_echo)
% - t0_echo - Starting time of the first gradient train (the first TE)
%   within a slice. The first gradient train within a TR will start at
%   t0_slice + t0_echo ;
% - GradRaster - The resolution of the waveform, all times are rounded to
%   integer multiples of GradRaster. On Siemens, GradRaster is 10 us.
%
% Outputs
% -------
% - G - The gradient waveform generated, with amplitude of 1, with each
%   point/sample being GradRaster long.
  
  NRamp = round(tRamp/GradRaster) ;
  NPlateau = round(tPlateau/GradRaster) ;
  NESP = round(tESP/GradRaster) ;
  NDeltaTE = round(DeltaTE/GradRaster) ;
  NTSlice = round(TSlice/GradRaster) ;
  NTR = round(TR/GradRaster) ;
  N0_slice = round(t0_slice/GradRaster) ;
  N0_echo = round(t0_echo/GradRaster) ;

  % Single gradient lobe
  % Note: Sampling at the edge of each gradient raster fits the model
  %       better than assuming the samples are at the center of each 
  %       gradient raster time.
  GRampUp = (0:(NRamp-1))/NRamp ;
  GRampDown = (NRamp:-1:1)/NRamp ;
  GPlateau = ones(1, NPlateau) ;
  NESPPad = NESP - (NRamp + NPlateau + NRamp) ;
  GLobe = [GRampUp, GPlateau, GRampDown, zeros(1, NESPPad)]  ;

  % Pair of gradient
  GPair = [GLobe, -GLobe] ;

  % Echo-train
  if (mod(ETL,2) == 0)
    GEchoTrain = repmat(GPair, 1, ETL/2) ;
  else
    GEchoTrain = [repmat(GPair, 1, (ETL-1)/2), GLobe] ;
  end    

  % Multi Echo
  if (NumTEs < 1)
    error('NumTEs must be at least 1') ;
  end
  if (NumTEs > 1)
    NTEPad = NDeltaTE - numel(GEchoTrain) ;
    GMultiTE = [zeros(1, N0_echo), ... (shift start)
                repmat([GEchoTrain, zeros(1, NTEPad)], 1, NumTEs-1), ...
                GEchoTrain] ;
  else
    GMultiTE = [zeros(1, N0_echo), ... (shift start)
                GEchoTrain] ;
  end

  % Multi slice
  if (NumSlices > 1)
    NSlicePad= NTSlice - numel(GMultiTE) ;
    GMultiSlice = [zeros(1, N0_slice), ... (shift start)
                   repmat([GMultiTE, zeros(1, NSlicePad)], 1, NumSlices-1), ...
                   GMultiTE] ;
  else
    GMultiSlice = [zeros(1, N0_slice), ... (shift start)
                   GMultiTE] ;
  end

  % Single TR
  NTRPad = NTR - numel(GMultiSlice) ;

  G = [GMultiSlice, zeros(1, NTRPad)] ;

end



