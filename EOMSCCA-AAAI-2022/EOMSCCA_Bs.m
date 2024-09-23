function label = EOMSCCA_Bs(Bs, nCluster)
% [1] Effcient One-Pass Multi-View Subspace 
% Clustering with Consensus Anchors-AAAI-2022
as = cellfun(@(x) size(x, 2), Bs);
anchor = min(as);
% anchor = nCluster;
d = (1)*nCluster ;
Bs = Bs(:);
nSmp = size(Bs{1}, 1);
Y = randi(nCluster, nSmp, 1);
[A,W,Z,iter,obj,alpha,label] = algo_qp(Bs, Y, d, anchor);
end