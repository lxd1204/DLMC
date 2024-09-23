function B = DiffusionBP_heat_n(B, varargin)
%**********************************************
%  [1] Diffusion Improves Graph Learning, NIPS, 2019
%  [2] Learning with Local and Global Consistency, NIPS, 2003
%  [3] The heat kernel as the pagerank of a graph-PNAS-2007
%**********************************************

[nSmp, nAnchor] = size(B);

param_names = {'eta', 'nNeighbor'};
param_default =  {5, 5};
[eid, errmsg, eta, nNeighbor] = getargs(param_names, param_default, varargin{:});
if ~isempty(eid)
    error(sprintf('DiffusionBP_heat:%s', eid), errmsg);
end

B = bsxfun(@times, B, 1./sqrt(max(sum(B, 1), eps)));
%*********************************************
% Step 1: The symmetric transition matrix
%*********************************************
G = B * B'; % (nm^2)
Ssym = (G + G')/2;
DSsym = 1./sqrt(max(sum(Ssym, 2), eps));
Gnorm = (DSsym * DSsym') .* Ssym;
Gnorm = (Gnorm + Gnorm')/2;
Gnorm = sparse(Gnorm);
%*********************************************
% Step 2: Heat Kernel diffusion with close-form
%*********************************************
L = eye(nSmp) - Gnorm;
L = full(L);
S = expm(- eta * L); % (m^3)
S = (S + S')/2;
% [V, D] = eig(S);
% S_half = V * diag(sqrt(diag(D)));
% BS = B * S_half; % (nm^2)
SB = S * B;

%*********************************************
% Step 3: Sparse it
%*********************************************
B = zeros(nSmp, nAnchor);
[~, Idx] = sort(SB, 2, 'descend');
Idx = Idx(:, 1:nNeighbor);
for iSmp = 1:nSmp
    idxa0 = Idx(iSmp, :);
    ad = SB(iSmp, idxa0);
    B(iSmp, idxa0) = EProjSimplex_new(ad);
end
B = bsxfun(@times, B, 1./sqrt(max(sum(B, 1), eps)));
end