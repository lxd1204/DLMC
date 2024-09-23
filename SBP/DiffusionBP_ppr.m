function B = DiffusionBP_ppr(B, varargin)
[nSmp, nAnchor] = size(B);

param_names = {'alpha', 'nNeighbor'};
param_default =  {20, 0.05};
[eid, errmsg, alpha, nNeighbor] = getargs(param_names, param_default, varargin{:});
if ~isempty(eid)
    error(sprintf('DiffusionBP_heat:%s', eid), errmsg);
end

B = bsxfun(@times, B, 1./sqrt(max(sum(B, 1), eps)));
mu = alpha/(1-alpha);
BTB = B' * B;
BTB = (BTB + BTB')/2;
IBTB =  (mu+1)/mu * eye(size(B, 2)) - BTB;
B2 = B * (IBTB \ BTB);
B2 = B + B2;
%*********************************************
% Step 3: Sparse it
%*********************************************
B = zeros(nSmp, nAnchor);
[~, Idx] = sort(B2, 2, 'descend');
Idx = Idx(:, 1:nNeighbor);
for iSmp = 1:nSmp
    idxa0 = Idx(iSmp, :);
    ad = B2(iSmp, idxa0);
    B(iSmp, idxa0) = EProjSimplex_new(ad);
end
% B = bsxfun(@times, B, 1./sqrt(max(sum(B, 1), eps)));
end