function [hijMat, mijMat, overlapHij, overlapMij, timeShiftSeconds] = ...
    connection_likelihood(x, sigmaH, maxH, sigmaM, delaySpread, relScaleDeviatationThd)
% Calculate H & M connection likelihoods
%
% Input:
% 1) x - struct array holding time & frequency measurements of whistle-traces
% 2) sigmaH - tuning parameter for the calculation of H.
% 3) maxH - maximal number of harmonics for the calculation of H.
% 4) sigmaM - tuning parameter for the calculation of M.
% 5) delaySpread - delay spread for multipath likelihood, in [s].
%
% Output:
% 1) hijMat - harmonic likelihoods, H.
% 2) mijMat - multipath likelihoods, M.
% 3) overlapHij - overlap for the calculation of H.
% 4) overlapMij - overlap for the calculation of M.
% 5) timeShiftSeconds - time shift between whistle-traces, in [s].
%

nSigDigits = 6; % Round times to nSigDigits
nTrackers = length(x);
debugFlag = 0;
tMax4Mij = 2*delaySpread;

% Initialize default values:
hijMat = zeros(nTrackers);
mijMat = zeros(nTrackers);

overlapHij = eye(nTrackers); % Default is 1 for the same tracker
overlapMij = eye(nTrackers); % Default is 1 for the same tracker
timeShiftSeconds = 0;
rationalRatios = calc_rational_ratio(maxH); % Ratinal frequency ratios for calculating H
dt = x(1).time(2)-x(1).time(1); % [s] time differene between trackers' samples
                 
for currentTrackNo = 1:nTrackers-1
    xi = x(currentTrackNo);
    t1 = xi.time(1);
    t2 = xi.time(end);
    for trackNo = currentTrackNo+1 : nTrackers
        xj = x(trackNo);
            if (t1 > xj.time(1)) && (t1 < xj.time(end)) || ...
                    (t1 <= xj.time(1)) &&  ( t2 >= xj.time(1))
                % Otherwise: used default values of Hij, Mij, etc.
                
                % --- Calculate harmonic likelihood, Hij ---
                [times4Correlation,ia,ib]  = intersect( round( xi.time, nSigDigits), ...
                    round(xj.time, nSigDigits)); % ia are the indices in the 'reference' tracker, and ib are the indeces in the 'test' tracker
                if numel(ia)>1 % verifay that overlap is greater than 1 (which produces 0 fitError)
                    overlapHij(currentTrackNo, trackNo) = numel(ia) / min( numel(xi.time), numel(xj.time));
                    %freqScale = mean(xi.freq(ia) ./xj.freq(ib)); % Original version - LS
                    freqScale = mean(xi.freq(ia)) / mean(xj.freq(ib));  % Modified version - to get symmetric Hij                            
                    scaleDeviation = min( abs(rationalRatios - freqScale)); % scaleDeviation is always positive
                    phiHij = scaleDeviation/freqScale < relScaleDeviatationThd; % phiH for nearly rational freqScale 
                    %phiHij = freqScale>=1/maxH && freqScale<=maxH; % phiH for any freqScale         
                    fitErrorHij = canberra_dist(xi.freq(ia), freqScale*xj.freq(ib));
                    hijMat(currentTrackNo, trackNo) = phiHij*exp(-(fitErrorHij/sigmaH)^2)*overlapHij(currentTrackNo, trackNo);
                    hijMat(trackNo, currentTrackNo) =   hijMat(currentTrackNo, trackNo); % Since hij=hji
                end
            end
            
            if abs(t1 - xj.time(1))<tMax4Mij ||  abs(t2 - xj.time(end))<tMax4Mij
                % --- Calculate multipath likelihood, Mij ---
                [fitErrorMij, relOverlap, timeShiftSeconds] = calc_time_shift_error(xi, xj, dt);
                phiM = abs(timeShiftSeconds) <= delaySpread;
                overlapMij(currentTrackNo, trackNo) = relOverlap;
                mijMat(currentTrackNo, trackNo) = phiM*exp(-(fitErrorMij/sigmaM)^2)*overlapMij(currentTrackNo, trackNo);
                mijMat(trackNo, currentTrackNo)  = mijMat(currentTrackNo, trackNo); % Since mij=mji
            end
            
            if debugFlag
                figure(104); clf
                subplot(2,1,1)
                plot(xi.time, xi.freq,'-o'); hold on
                plot(xj.time,  xj.freq ,'-o')
                plot(xj.time,  freqScale*xj.freq ,'-+')
                legend(['i=' num2str(currentTrackNo)], ['j=' num2str(trackNo)], 'reconstruction')
                title(['H_{ij}=' num2str( hijMat(currentTrackNo, trackNo) ) ',F_s=' num2str(freqScale) ...
                    ', \epsilon_H=' num2str( fitErrorHij) ', O_H=' num2str(overlapHij(currentTrackNo, trackNo))])
                subplot(2,1,2);
                plot(xi.time, xi.freq,'-o'); hold on
                plot(xj.time,  xj.freq ,'-o')
                plot(xj.time + timeShiftSeconds,  xj.freq ,'-+')         
                title(['M_{ij}=' num2str( mijMat(currentTrackNo, trackNo) ) ',\DeltaT=' num2str(timeShiftSeconds) ...
                    '[s], \epsilon_M=' num2str(fitErrorMij) ', O_M=' num2str(overlapMij(currentTrackNo, trackNo))])
            end            
    end
end



function fitError = canberra_dist(q1, q2)
% Calculate normalized Canberra distance between q1 and q2:
fitError = sum( abs(q1 - q2 ) ./ ( abs(q1) + abs(q2) ) );
fitError = fitError/numel(q1); % Normalize by the length of the compared sequences


function  [fitErrorMij, relOverlap, timeShiftSeconds] = calc_time_shift_error(xi, xj, dt)
if length(xi.time)>length(xj.time)
    x1 = xi; x2 = xj;
    timeShiftSign = 1;
else
    x1 = xj; x2 = xi;
    timeShiftSign = -1;
end
f1ZeroPadded = [zeros(1, length(x2.freq)-1)   (x1.freq(:))'    zeros(1, length(x2.freq)-1)];
%f1ZeroPadded = [zeros(1, length(x2.freq)-1)   x1.freq    zeros(1, length(x2.freq)-1)]; % Original line
f2ZeroPadded = zeros(size(f1ZeroPadded)); f2ZeroPadded(length(x2.freq):2*length(x2.freq)-1) = x2.freq;
timeShiftVec = -(length(x2.freq)-1):(length(x1.freq)-1);
fitErrorVec = ones(size(timeShiftVec)) * 1e6;
nValidIndVec = zeros(size(timeShiftVec));
for timeShiftInd = 1:length(timeShiftVec)
    f2ZpShifted = circshift(f2ZeroPadded, timeShiftVec(timeShiftInd));
    validInd =  f2ZpShifted>0 & f1ZeroPadded>0;
    fitErrorVec(timeShiftInd) = canberra_dist(f1ZeroPadded(validInd), f2ZpShifted(validInd));
    nValidIndVec(timeShiftInd) = nnz(validInd) /min( numel(x1.freq), numel(x2.freq));
    %nValidIndVec(timeShiftInd) = (nnz(validInd)>1)*nnz(validInd) /min( numel(xi.freq), numel(xj.freq)); % verifay that overlap is greater than 1 (which produces 0 fitError)
end
[~, bestTimeShiftInd] = min(fitErrorVec./nValidIndVec);
fitErrorMij = fitErrorVec(bestTimeShiftInd);
timeShiftSeconds = timeShiftSign*timeShiftVec( bestTimeShiftInd )*dt + xi.time(1) - xj.time(1);
relOverlap = nValidIndVec(bestTimeShiftInd);