function LLR1 = dilute_LLR(LLR)
% Dilute LLR matrix to only maximal likelihoods
LLR1 = zeros( size( LLR ) );
for ii = 1:size(LLR1,1)
    [valOfMaxLLR, indOfMaxLLR] = max(LLR(ii,:));
    LLR1(ii, indOfMaxLLR) = valOfMaxLLR;
end
LLR1 = 0.5*(LLR1 + LLR1');