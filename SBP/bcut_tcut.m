function [U_normalized, V_normalized] = bcut_tcut(B, nCluster)
% 
% Reference
% [2] Segmentation Using Superpixels: A Bipartite Graph Partitioning Approach
% 
[nSmp, nAnchor] = size(B);

% build anchor-anchor graph
dx = sum(B,2);
dx(dx==0) = eps; 
Dx = sparse(1:nSmp,1:nSmp,1./dx); %nSmp * nSmp
clear dx
Wy = B'*Dx*B;  % nAnchor * nAnchor

% normalized affinity matrix
d = sum(Wy,2);
d(d==0) = eps;
D = sparse(1:nAnchor,1:nAnchor,1./sqrt(d)); % nAnchor * nAnchor
clear d
nWy = D*Wy*D; 
clear Wy
nWy = (nWy+nWy')/2;  % nAnchor * nAnchor

%%% compute Ncut eigenvectors
% computer eigenvectors
[evec,eval] = eig(full(nWy)); 
clear nWy % use eigs for large superpixel graphs  
[~, idx] = sort(diag(eval), 'descend');
Ncut_evec = D*evec(:,idx(1:nCluster)); % nAnchor * nCluster
clear D

% compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
evec = Dx * B * Ncut_evec; 
if nargout > 1
    V_normalized = bsxfun( @rdivide, Ncut_evec, sqrt(sum(Ncut_evec.*Ncut_evec,2)) + 1e-10 );
end

clear B Dx Ncut_evec

% normalize each row to unit norm
U_normalized = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );
end