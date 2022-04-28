function [kappaij, LLR, SNR] = detect_continuity(y, fs, fStart, fEnd, sigmaKappa, nTests, ...
    minSearchRange, plotFlag, rngSeed)
% Detect continuity of a chirp componenet in the signal y, between frequencies [fStart  fEnd]
%
% Input:
% 1) y - input signal.
% 2) fs - sampling frequency, in [Hz].
% 3) fStart - chirp frequency at the start of the signal, in [Hz].
% 4) fEnd - chirp's frequency at the end of the signal, in [Hz].
% 5) sigmaKappa - tunning parameter for the calculatin of kappaij
% 6) nTests - number of tests for the LLR.
% 7) minSearchRange -  minimal search range around the asuumed location of the peak in the fractional domain (in samples).
% 8) plotFlag - 1/0 - do/don't make plots
% 9) rngSeed - seed for random numbers generator, to enable reproduction of
%       scenarios. Leave empty if you don't want to use a pre-determined seed.
%
% Output:
% 1) Cij - continuity likelihood
% 2) LLR - log likelihood ratio of the continuity test
% 3) SNR - estimated SNR of detected peak
%
%  
% detect_continuity
%               |
%               |---- fracF
%               |
%               |---- get_peak_range
%               |
%               |---- get_params
%               |
%               |---- calc_kappaij
%

if ~isempty(rngSeed)
    rng('default');
    rng(rngSeed);
end
if mod(length(y), 2) % Make sure that length of signal is even
    y = y(1:end-1);
end

chirpRate = (fEnd - fStart) / (length(y)/fs); % [Hz]/[s]
optimalAlpha =  -(2/pi)*atan((fs^2/length(y))./chirpRate);

% -- Calculate FrFT for the slope corresponding to chirpRate ----
Y = fracF(y, optimalAlpha);
[expectedPeakLoc, searchRange] = get_peak_range(length(Y), [fStart  fEnd], optimalAlpha, fs, minSearchRange);
ind1 = max(round( expectedPeakLoc - searchRange), 1);
ind2 = min(round( expectedPeakLoc + searchRange), length(Y));
[maxPeak, peakLoc] = max(abs(Y(ind1:ind2)));
peakLoc = peakLoc + ind1 - 1;

% --- Calculate SNR of detected peak
 [~, SpectrogramParams] = get_params;
 medfiltBandwidth = SpectrogramParams.medfiltLen * fs/SpectrogramParams.nfft; % [Hz]
 df = fs/length(Y); 
 medfiltLen = round(medfiltBandwidth/df);
 ind1 = max( round(peakLoc - medfiltLen/2), 1);
 ind2 = min( round( peakLoc  + medfiltLen/2), length(Y) );
 SNR = 20*log10(maxPeak) - 20*log10(median(abs(Y(ind1:ind2))));
 
% --- Make plot 
if plotFlag
    figure
    semilogy(abs(Y))
    hold on
    plot(peakLoc,maxPeak,'xr');
    plot(round(expectedPeakLoc), abs(Y(round(expectedPeakLoc))), 'or');
    legend('\alpha_{opt}', 'detected peak','expected peak');
    title(['SNR = ' num2str(SNR) 'dB'])
end


% -- Calculate FrFT for nTests random slopes ----
alphaVec = rand(nTests, 1)*2-1;
maxPeak1 = zeros(nTests, 1);
peakLoc1 = zeros(nTests, 1);
for ii = 1:length(alphaVec)
    Y1 = fracF(y, alphaVec(ii));
    [maxPeak1(ii), peakLoc1(ii)] = max(abs(Y1));
    if plotFlag && ii<=5
        legendCell{ii+3} = ['\alpha = ' num2str(alphaVec(ii))];
        semilogy(abs(Y1))
    end
end
if plotFlag
    legend(legendCell)
end

% -- Log-likelihood ratio test ---
LLR = log(maxPeak) - sum( log( maxPeak1) )/nTests;
kappaij = calc_kappaij(LLR, sigmaKappa/sqrt(nTests)); 
