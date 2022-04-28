function kappaij = calc_kappaij(LLR, sigmaKappa)
% Calculate continuity likelihood
kappaij = exp(LLR/sigmaKappa^2) ./ (1 + exp(LLR/sigmaKappa^2));
