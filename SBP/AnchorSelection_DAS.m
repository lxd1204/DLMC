function [anchors, ind_all, score_all] = AnchorSelection_DAS(X, nAnchor)
% 
% [1] Multi-view Clustering A Scalable and Parameter-free Bipartite Graph Fusion Method, TPAMI, 2022
% 
% 
% Code is provided by Author Han Zhang
% 
[nSmp, nFea] = size(X);

% *******************************************************************************************
% Transformation 
% 
% where Tra() denotes the operation that transforms the raw features to be non-negative. 
% More specifically, when features are not non-negative, we translate the features in each
% dimension by subtracting the minimum value within each dimension. 
% The motivation of this operation is that the features of samples in the same cluster generally 
% have smaller differences than features of samples in difference clusters,
% and samples in the same cluster naturally obtain similar scores computed by summing the feature values.
% 
% *******************************************************************************************
vm = min(X, [], 1);
X = bsxfun(@minus, X, vm);
score = sum(X, 2);

score_all = zeros(nSmp, nAnchor);
ind_all = zeros(nAnchor, 1);
score_all(:, 1) = score/max(score);
[~, ind_all(1)] = max(score_all(:, 1));
for iAnchor = 2:nAnchor
   score_all(:, iAnchor) = score_all(:, iAnchor - 1) .* (ones(nSmp, 1) - score_all(:, iAnchor -1));
   score_all(:, iAnchor) = score_all(:, iAnchor)/max(score_all(:, iAnchor));
   [~, ind_all(iAnchor)] = max(score_all(:, iAnchor));
end
ind_all = sort(ind_all, 'ascend');
anchors = X(ind_all, :);
end