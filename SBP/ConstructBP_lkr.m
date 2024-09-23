function B = ConstructBP_lkr(X, anchors, varargin)
% Input
%         X: nSmp * nFea
%         anchors: nAnchor * nFea
%         nNeighbor: row sparsity
% Output
%         B: nSmp * nAnchor
%
%

[nSmp, nFea] = size(X);
[nAnchor, nFea] = size(anchors);

param_names = {'nNeighbor'};
param_default =  {5};
[eid, errmsg, nNeighbor] = getargs(param_names, param_default, varargin{:});
if ~isempty(eid)
    error(sprintf('ConstructBP_lkr:%s', eid), errmsg);
end

aa = sum(X.^2, 2);
bb = sum(anchors.^2, 2);
Kb = X * anchors';
D = bsxfun(@plus, aa, bb') - 2 * Kb;
D (D < 1e-10) = 1e10;  % anchor_j = x_i
[D2, Idx] = sort(D, 2); % sort each row
Idx = Idx(:, 1:nNeighbor);

Ka = anchors * anchors';
Ka = (Ka + Ka')/2;
Ik = eye(nNeighbor);
e = ones(1, nNeighbor);
z = zeros(nNeighbor, 1);
e2 = ones(nNeighbor, 1);
options = [];
options.Display = 'off';
B = zeros(nSmp, nAnchor);
for iSmp = 1:nSmp
    idx = Idx(iSmp, :);
    ki = Kb(iSmp, idx)';
    Kii = Ka(idx, idx') + Ik;
    v = quadprog(Kii, -ki, [], [], e, 1, z, e2, [], options);
    B(iSmp, idx) = v';
end

end