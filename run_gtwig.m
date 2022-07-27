% Run GTWIG example
%
%
% *IMPORTANT NOTE:
% This codes needs an external implementation of fractional Fourier transform. I used the 
% function fracF.m that can be downloaded from:
% http://www.ee.bilkent.edu.tr/~haldun/wileybook.html
% (Matlab code for fast computation of the fractional Fourier transform)
% Alternately, you can use other implementations, or disable the trace-continuity calculation by
% setting:  AlgParams.calcContinuityFlag=0


%% --- Data files ---
audioFileName = '2021_06_30_18_33_31_1850_short.flac' % 48 Original Traces --> 39 Segmented Traces
%audioFileName = '2021_06_26_18_22_24_700.flac'; % 165 original traces --> 161 segmented traces (for AlgParams.calcContinuityFlag=1)
trackersFileName = [audioFileName(1:end-5) '_GT.mat'];
FileParams = struct(...
        'audioPath',  '.\DEMO_DATA', ... % Directory of wav files
    'trackersPath', '.\DEMO_DATA'); % Directory for whistle-traces

%% --- Parameters for the algorithm ---
AlgParams = struct( ...
    'delaySpread', 0.27, ... % [s]  Delay to consider multipath
    'calcContinuityFlag', 0, ... % *see note above
    'maxDeltaT', 0.1, ... % [s] Maximal time gap between traces for continuity test
    'maxH', 5, ... % Maximal number of harmonics assumed
    'minSearchRange', 5, ...  % minimal search range around the asuumed location of the peak in the fractional domain (in samples).
    'nTests4LLR', 30, ... % Number of tests for the log-likelihood ratio
    'relScaleDeviatationThd', 0.01, ... % Relative frequency deviation to consider two frequencies ahrmonics or not
    'sigmaH', 0.0003, ... % Value used in paper is 0.0003
    'sigmaM', 0.00015); % Value used in paper is 0.0001
    
  
%% -- run G-TWIG ---
[x, idxClustered] =  gtwig(audioFileName, trackersFileName, AlgParams, FileParams);


%% -- Plot spectrogram of audio data ----
[y, fs] = audioread([FileParams.audioPath filesep audioFileName]);
[s,f,t] = spectrogram(y(:,1), 4096, round(4096*0.8), 4096, fs); 
figHandle = figure; imagesc(t, f/1000, 20*log10(abs(s))); colormap gray
xlabel('time [s]'); ylabel('frequency [kHz]'); axis xy

%% -- Plot clusters ---
plot_clusters(x, idxClustered, figHandle, 1);
