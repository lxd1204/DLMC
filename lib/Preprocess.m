function X = Preprocess(X, varargin)
[nSmp, nFea] = size(X);

param_names = {'type', 'norm'};
param_default =  {1, 2};
[eid, errmsg, type, norm] = getargs(param_names, param_default, varargin{:});
if ~isempty(eid)
    error(sprintf('Preprocess:%s', eid), errmsg);
end

switch type
    case 1
        xNorm = max(1e-14, full(sum(abs(X).^norm, 2) ) );
        X = spdiags(xNorm.^-(1./norm), 0, nSmp, nSmp) * X;
    case 2
end
end