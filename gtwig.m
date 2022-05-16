function [x, idxClustered] = gtwig(wavFileName, trackersFileName, AlgParams, FileParams)
% Cluster whistle traces by the  G-TWIG algorithm, presented in:
% D. Kipnis and R. Diamant, "Graph-Based Clustering of Dolphin Whistles." IEEE/ACM
% Transactions on Audio, Speech, and Language Processing 29 (2021): 2216-2227.
%
% Code written by Dror Kipnis, 2020-2021
%
%
% *IMPORTANT NOTE:
% This codes needs an external implementation of fractional Fourier transform. I used the
% function fracF.m that can be downloaded from:
% http://www.ee.bilkent.edu.tr/~haldun/wileybook.html
% (Matlab code for fast computation of the fractional Fourier transform)
% Alternately, you can use other implementations, or disable the trace-continuity calculation by
% setting:  AlgParams.calcContinuityFlag=0
%
%
% Input:
% 1) wavFileName - name of audio file.
% 2) trackersFileName - name of file that contain whistle-traces.
% 3) AlgParams - parameters of the algorithm.
% 4) FileParams - contain directoreis for data files (wavFileName & trackersFileName).
%
% Output:
% 1) x - whistle traces loaded from data file trackersFileName.
% 2) idxClustered - cluster number associated with each whistle-trace in x.
%
%
%
% Functions tree:
%
%      gtwig
%          |
%          |---- get_params
%          |
%          |---- idx2Ck
%          |
%          |---- connection_likelihood
%          |                |
%          |                |---- calc_rational_ratio
%          |
%          |---- detect_continuity
%          |                |
%          |                |---- fracF*
%          |                |
%          |                |---- get_peak_range
%          |                |
%          |                |---- get_params
%          |                |
%          |                |---- calc_kappaij
%          |
%          |---- calc_kappaij
%          |
%          |---- dilute_LLR
%          |
%          |---- graph_spectral_clustering
%          |
%          |--- Ck2idx
%          |
%          |---- plot_clusters
%                       |
%                       |---- get_trackers_no


% -- Load trackers ---
TrackersData = load([FileParams.trackersPath filesep trackersFileName]);
x = TrackersData.GT;
clear TrackersData;
nTrackers = length(x);

if nTrackers>1
    % -- Load audio ---
    if AlgParams.calcContinuityFlag % Audio file is required only for the continuity calculation
    	disp(['Loading file: ' FileParams.audioPath  filesep wavFileName]);
    	[sig, fs] = audioread([FileParams.audioPath  filesep wavFileName]);
    end

    %% -- Calculate H & M -------
    disp('Calculating Hij & Mij');
    [Hij, Mij] = connection_likelihood(x, AlgParams.sigmaH, AlgParams.maxH, AlgParams.sigmaM, ...
        AlgParams.delaySpread, AlgParams.relScaleDeviatationThd);


    %% -- Calculate Kappaij ---------
    % Identify candidate traces, take signal section between their ends, and
    % perform continuity detection
    if AlgParams.calcContinuityFlag
        sigmaKappa = 1/sqrt(AlgParams.nTests4LLR);
        LLRThd = 9/ AlgParams.nTests4LLR;
        rngSeed = []; % 0
        LLR = zeros(nTrackers, nTrackers);
        disp('Calculating LLR');
        for ii = 1:nTrackers
            tQi = x(ii).time(end);
            for jj = 1:nTrackers
                t1j = x(jj).time(1);
                if tQi<t1j && (t1j-tQi) <= AlgParams.maxDeltaT
                    disp(['Testing trackers ' num2str(ii) ', ' num2str(jj)])
                    y = sig(floor(tQi*fs):ceil(t1j*fs));
                    fStart = x(ii).freq(end);
                    fEnd = x(jj).freq(1);
                    [~,  LLR(ii, jj) ] = detect_continuity(y, fs, fStart, fEnd, LLRThd, AlgParams.nTests4LLR, ...
                        AlgParams.minSearchRange, 0, rngSeed);
                    LLR(jj, ii) = LLR(ii, jj);
                end
            end
        end
        %save(continuityFileName, 'LLR'); % Calculation takes time, so you can save the LLR and reload it later

        kappaij = calc_kappaij(dilute_LLR(LLR) - LLRThd, sigmaKappa);
        kappaij(LLR <=0) = 0;
    else
        kappaij = zeros(nTrackers, nTrackers);
    end

    %% -- Calculate W_ij ------
    wijMat = max(max(Hij, Mij), kappaij); % ML


    %% --- Clustering ------
    D = diag( sum(wijMat,2)); % Option 1: h_ii==0
    Ck =  graph_spectral_clustering(D-wijMat); % Unnormalized graph Laplacian
    idxClustered = Ck2idx(Ck);
else
    if nTrackers==1
        idxClustered = 1;
    else % nTrackers==0
        idxClustered = [];
    end
end
