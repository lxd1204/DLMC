function paramCell = EMKMC_build_param(anchor_sizes, gammas)

if ~exist('anchor_sizes', 'var')
    anchor_sizes = [[50, 100]; [50, 100]];
end

if ~exist('gammas', 'var')
    gammas = [1.6];
end

nParam = size(anchor_sizes, 1) * length(gammas);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:size(anchor_sizes, 1)
    for i2 = 1:length(gammas)
        param = [];
        param.anchors = anchor_sizes(i1, :);
        param.gamma = gammas(i2);
        idx = idx + 1;
        paramCell{idx,1} = param;
    end
end
end