function label = FastMICE_Bs(Bs, nCluster)
nView = length(Bs);
nSmp = size(Bs{1}, 1);
idx = zeros(nSmp, nView);
for iView = 1:nView
    idx(:, iView) = Tcut_for_bipartite_graph(Bs{iView}, nCluster);
end
label = FastMICE_ConsensusFunction(idx, nCluster, 100, 10);
end