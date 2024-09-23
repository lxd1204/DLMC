function [anchors, ind_all, score_all] = AnchorSelection_ALG(A, nAnchor)
% 
% [1] Large-Scale Clustering With Structured Optimal Bipartite Graph
% 
% 
% Code is provided by Author Han Zhang
% 
[nSmp, nSmp2] = size(A);

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

score_all = zeros(nAnchor, nSmp);
ind_all = zeros(nAnchor, 1);
d = sum(A, 1);
score_all(1, :) = d;
[~, ind_all(1)] = max(score_all(:, 1));
for iAnchor = 2:nAnchor
   score_all(iAnchor, :) = score_all(iAnchor - 1, :) .* (ones(1, nSmp) - score_all(iAnchor -1, :)) .* d;
   score_all(iAnchor, :) = score_all(iAnchor, :)/max(score_all(iAnchor, :));
   [~, ind_all(iAnchor)] = max(score_all(iAnchor, :));
end
ind_all = sort(ind_all, 'ascend');
anchors = A(:, ind_all);
end