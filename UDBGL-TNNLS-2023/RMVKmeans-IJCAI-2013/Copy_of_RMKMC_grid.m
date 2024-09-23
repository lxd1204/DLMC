function [resAIOCell] = RMKMC_grid(Xlow, Xup, Y, gammaCandi, nRepeat, mRepeat)
[nSmp, nFea] = size(Xlow);
nCluster = length(unique(Y));

% nRepeat = 10;
% mRepeat = 1;
maxIter = 100;



paramCell = RMKMC_build_param(gammaCandi, maxIter);
nParam = length(paramCell);
resAIOCell = cell(1, nParam);

inXCell = cell(2,1);
inXCell{1,1} = Xlow';
inXCell{2,1} = Xup';
for iParam = 1:nParam
    param = paramCell{iParam};
    resAIO = zeros(nRepeat, 6);
    
    for iRepeat = 1:nRepeat
        disp(['RMKMC: ', num2str(iParam), ' out of ', num2str(nParam), ';', num2str(iRepeat), ' out of ', num2str(nRepeat)]);
        tStart = tic;
        inPara = [];
        inPara.maxIter = maxIter;
        inPara.thresh = 1e-5;
        inPara.numCluster = nCluster;
        inPara.r = param.gamma;
        label = randi(nCluster, nSmp, 1);
        inG0 = sparse([1:nSmp]', label, 1, nSmp, nCluster, nSmp);
        [ outG0, ~, ~, ~, objHistory, ~ ] = weighted_robust_multi_kmeans( inXCell, inPara, inG0 );
        [~, label] = find(outG0);
        
        [ACC, MIhat, Purity, ARI] = ClusteringMeasure(Y, label);
        tElapsed = toc(tStart);
        objEnd = objHistory(end);
        resAIO(iRepeat, :) = [ACC, MIhat, Purity, ARI, tElapsed, objEnd];
    end
    
    resAIOCell{iParam} = resAIO;
end