function [resAIOCell] = RMKMC_grid(inXCell, Y, gammaCandi, nRepeat, mRepeat)
[nSmp, nFea] = size(inXCell{1});
nCluster = length(unique(Y));

% nRepeat = 10;
% mRepeat = 1;
maxIter1 = 100;
maxIter2 = 1;



paramCell = RMKMC_build_param(gammaCandi, maxIter1);
nParam = length(paramCell);
resAIOCell = cell(1, nParam);


for iParam = 1:nParam
    param = paramCell{iParam};
    resAIO = zeros(nRepeat, 6);
    
    for iRepeat = 1:nRepeat
        disp(['RMKMC: ', num2str(iParam), ' out of ', num2str(nParam), ';', num2str(iRepeat), ' out of ', num2str(nRepeat)]);
        tStart = tic;
        [label, ~, ~, objHistory] = RMKMC(inXCell, param.gamma, nCluster, maxIter1, maxIter2, mRepeat);
        [ACC, MIhat, Purity, ARI] = ClusteringMeasure(Y, label);
        tElapsed = toc(tStart);
        objEnd = objHistory(end);
        resAIO(iRepeat, :) = [ACC, MIhat, Purity, ARI, tElapsed, objEnd];
    end
    
    resAIOCell{iParam} = resAIO;
end