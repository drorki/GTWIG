function [expectedPeakLoc, searchRange] = get_peak_range(sigLen, fRange, a, ...
    fs, minSearchRange)
% Return the estimated location of the peak on the fractional-grid, and the desired
% search range around it.
% The estimated location of the peak is based on [Capus2000], [Capus2003] and [Cowell2010]
% The search range is calculated such that the uncertainty area is conserved.
%
% Input:
% 1) sigLen - total number of samples.
% 2) fRange - frequency range: [fStart  fEnd], in [Hz].
% 3) a - order of the fractional Fourier transform (FrFT).
% 4) fs - sampling frequency, in [Hz].
% 5) minSearchRange - minimal search range, in samples.
% 
% Output:
% 1) expectedPeakLoc - expected location of peak on the fractional frequency axis, in samples.
% 2) searchRange - uncertainty range around the peak, in samples.

Nf = sigLen*abs(sin(a*pi/2)); % Based on [Capus2000]
f0 = fRange(1);
deltaF = diff(fRange); % fEnd - fStart
mAlpha = f0*sigLen*sin( a*pi/2)/fs;  % Offset from center of fractional domain -  based on [Capus2003]
mCenter = deltaF*sigLen*0.5*sin( a*pi/2)/fs; % Offset between fractional and frequency axes - [Cowell2010]
expectedPeakLoc = Nf/2 + mAlpha  + mCenter; % Good option without fftshift
searchRange = abs( sigLen*deltaF*sin(a*pi/2)/fs );
searchRange = max(searchRange, minSearchRange);