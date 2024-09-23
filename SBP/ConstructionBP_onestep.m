function [Z, A] = ConstructionBP_onestep(X, m, beta)
% 
% Eq3 in [1]
% 
% [1]Align then Fusion: Generalized Large-scale Multi-view 
% Clustering with Anchor Matching Correspondences, NIPS, 2022
% 
% Input
%     X, nSmp * nFea
%     m, the number of anchors
%     beta, the regularization, a large beta leads to more sparse Z
% Output
%     Z, nSmp * m, bipartite graph
%     A, nFea * m, anchors
%

if ~exist('beta', 'var')
    beta = 0.1;
end
maxIter = 50 ; % the number of iterations
[nSmp, nDim] = size(X);

rand('twister',5489);
[~, A] = litekmeans(X, m, 'MaxIter', 100, 'Replicates', 10);
A = A';
X = mapstd(X',0,1); % turn into d*n

AX = (A'*X);
idx = 1:m;
Z = zeros(m, nSmp);
for iSmp = 1:nSmp
    Z(idx,iSmp) = EProjSimplex_new(AX(idx,iSmp));
end

iter = 0;
objHistory = [];
converges = false;
while ~converges
    iter = iter + 1;
    %*****************************
    % Optimize Anchors
    %*****************************
    if size(A,1) < size(A,2)
        A = X*Z'*pinv(Z*Z');
    else
        C = X*Z';
        [U,~,V] = svd(C,'econ');
        A = U*V';
    end
    
    %*****************************
    % Optimize Bipartite Graph
    %*****************************
    AX = (A'*X)/(1+beta);
    Z = zeros(m, nSmp);
    for iSmp = 1:nSmp
        Z(idx,iSmp) = EProjSimplex_new(AX(idx,iSmp));
    end
    
    o1 = norm(X - A * Z,'fro')^2;
    o2 = norm(Z,'fro')^2;
    obj = o1 + beta * o2;
    objHistory = [objHistory; obj];%#ok
    
    if (iter > 10 && abs(objHistory(end-1)- objHistory(end))/objHistory(end-1) < 1e-3) || iter > maxIter
        converges = true;
    end
end
Z = Z';
end