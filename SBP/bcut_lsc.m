function U_normalized = bcut_lsc(Z, nCluster)
% Z: nSmp * mLandmark
% Reference:
% mySVD.m
% LSC_eigen.m
% [1] Large Scale Spectral Clustering Via Landmark-Based Sparse Representation
%

[nSmp, nAnchor] = size(Z);

%*********************************************
% Row normalization
%*********************************************
Z = bsxfun(@times, Z, max(sum(Z, 2), 1e-2));

%*********************************************
% Column normalization
%*********************************************
feaSum = full(sqrt(sum(Z, 1)));
feaSum = max(feaSum, 1e-12);
Z = Z./feaSum(ones(size(Z,1),1),:);
U = mySVD(Z, nCluster + 1);
U(:, 1) = [];
U_normalized = U./repmat(sqrt(sum(U.^2,2)),1,size(U, 2));
end