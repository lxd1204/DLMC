function cvi_score = ComputeICVIs(X, label, cvis)
[nSmp, nFea] = size(X);
nSmp2 = length(label);
if ~exist('cvis', 'var')
    cvis = {'CalinskiHarabasz', 'DaviesBouldin', 'silhouette'};
end
nCVIs = length(cvis);
cvi_score = zeros(nCVIs, 1);
for iCVI = 1:nCVIs
    cvi = cvis{iCVI};
    eva = evalclusters(X, label, cvi);
    cvi_score(iCVI) = eva.CriterionValues;
end
end