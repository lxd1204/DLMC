function [Xa] = AnchorSelection_lkm(X, varargin)
%This function is to select anchors
% 
[nSmp, nFea] = size(X);

param_names = {'nAnchor', 'maxIter', 'nRepeat'};
param_default =  {1, 100, 10};
[eid, errmsg, nAnchor, maxIter, nRepeat] = getargs(param_names, param_default, varargin{:});
if ~isempty(eid)
    error(sprintf('AnchorSelection_lkm:%s', eid), errmsg);
end

if nAnchor > nSmp
    nAnchor = nSmp;
    Xa = X;
else
    [~, Xa] = litekmeans(X, nAnchor, 'MaxIter', maxIter, 'Replicates', nRepeat);
end

end