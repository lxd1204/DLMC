function label = LMVSC_Bs(Bs, nCluster)
G = cell2mat(Bs);
[U,~,~] = mySVD(G, nCluster); 

rand('twister',5489);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,size(U,2));
label = litekmeans(U_normalized, nCluster, 'MaxIter', 100, 'Replicates', 10);
end