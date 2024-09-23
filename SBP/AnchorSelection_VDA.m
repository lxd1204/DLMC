function [anchor, ind_all, score_all] = AnchorSelection_VDA(X, nAnchor)
% 
% [1] Tensorized Bipartite Graph Learning for Multi-View Clustering
% 
% 
% Code is provided by Author 
% 

[nSmp, nFea] = size(X);
X_std = std(X, [], 2);

score_all = zeros(nSmp, nAnchor);
ind_all = zeros(nAnchor, 1);

score = X_std.^2;
score_all(:, 1) = score/max(score);
[~, ind_all(1)] = max(score_all(:, 1));
for iAnchor = 2:nAnchor
    for j=1:nSmp
        A_1 = score_all(ind_all(iAnchor -1), iAnchor-1);
        A_2 = score_all(j, iAnchor-1);
        Co(j, :) = (1 + norm(A_1-A_2,2)^(0.5))^(-1);
    end
    pho = Co/max(Co);
    score_all(:, iAnchor) = score_all(:, iAnchor-1).*(ones(nSmp,1)-pho);
    score_all(:, iAnchor) = score_all(:, iAnchor)/max(score_all(:, iAnchor));
    [~, ind_all(iAnchor)] = max(score_all(:, iAnchor));
end
ind_all = sort(ind_all, 'ascend');
anchor = X(ind_all, :);
end