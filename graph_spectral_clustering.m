function Ck = graph_spectral_clustering(L)
% Perform spectral clustering using graph Laplacian matrix L=D-W
% Stages are:
% a) Find eigenvectors and eigenvalues of the graph Laplacian L=D-W
% b) Identify eigenvectors corresponding to the K 2nd smallest eigenvalues, 
%     that contain all the nodes.
% c) Cluster by assigning each item to the cluster on which it has the largest projection
%
% Input:
% L - uunormalized graph laplacian matrix,  L=D-W, where W is the weighted
%   adjacency matrix, and D is a diagonal matrix whose (i,i) entry is the sum of 
%   similarities of the i-th whistle-trace to all other whistle-traces
%
% Output:
% Ck - Clustring solution, where each column corresponds to a cluster.
%

[U, S, V] = svd(L); % U = left singular vecotrs in clumns, S = singular values in diagonal, V = right singular vectors in clumns

% -- Evaluate model order K - identify the fiert eigenvector that contain each item ---
nItems = size(L, 1);
firstInd = ones(nItems, 1);
tol = 1/nItems;
[~, ind] = sort(diag(S),'ascend'); % eigenvalues sorted in increasing order
V1 = U(:, ind); % S sorted in increasing ordr of eigenvalues
for jj = 1:nItems
    indOfAssumedCluster = find(abs( V1(jj, 2:end) ) > tol, 1 );
    if ~isempty(indOfAssumedCluster)
        firstInd(jj) = indOfAssumedCluster;
    end
end
K = max(firstInd);

Ck = zeros(size(V1)); % Initialize according to number or eigenvectors
% --- form clusters by MAP -----
ind2Keep = zeros(1, nItems); % Indicate cluster into whish trackers are assigned
for jj = 1:nItems
    [~, k] = max( abs(V1(jj,2:(2+K-1)) ) ); 
    ind2Keep(jj) = k+1;
    Ck(jj,  ind2Keep(jj)) = 1;
end
% --- Prune empty clusters ----
 Ck = Ck(:, unique(ind2Keep));
 