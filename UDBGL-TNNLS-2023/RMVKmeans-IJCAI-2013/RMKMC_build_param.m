function paramCell = RMKMC_build_param(gammaCandi, maxIterCandi)

if ~exist('gammaCandi', 'var') || isempty(gammaCandi)
    gammaCandi = [1];
end

if ~exist('maxIterCandi', 'var') || isempty(maxIterCandi)
    maxIterCandi = 100;
end

nParam = length(gammaCandi) * length(maxIterCandi);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(gammaCandi)
    for i2 = 1:length(maxIterCandi)
        param = [];
        param.gamma = gammaCandi(i1);
        param.maxIter = maxIterCandi(i2);
        idx = idx + 1;
        paramCell{idx,1} = param;
    end
end