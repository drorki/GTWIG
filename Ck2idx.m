function idx = Ck2idx(Ck)
% Return cluster indices for the clustring assignment matrix Ck.
% 
% Input:
% 1) Ck - NxK binary matrix containing indicators for which of the N items
% belong to the k-th cluster, k=1,..K

N = size(Ck, 1); % Number of items
idx = zeros(N, 1);
for ii = 1:N
    idx(ii) =  find(Ck(ii,:));
end
