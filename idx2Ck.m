function Ck = idx2Ck(idx)
% Return matrix with cluster vectors in is columns, for the clustering solution idx
% 
% Input:
% 1) idx - cluster index corresponding to each item
%
% Output:
% 1) Ck - NxK binary matrix containing indicators for which of the N items
% belong to the k-th cluster, k=1,..K

N = length(idx);
K = numel(unique(idx));
Ck = zeros(N, K);
for itemNo = 1:N
    Ck(itemNo, idx(itemNo)) = 1;
end
