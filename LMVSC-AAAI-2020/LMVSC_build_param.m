function paramCell = LMVSC_build_param(anchor_sizes, alphas)

if ~exist('anchor_sizes', 'var')
    anchor_sizes = [50, 100];
end

if ~exist('alphas', 'var')
    alphas = [0.001, 0.01, 0.1, 1, 10];
end

nParam = length(anchor_sizes) * length(alphas);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(anchor_sizes)
    for i2 = 1:length(alphas)
        param = [];
        param.nAnchor = anchor_sizes(i1);
        param.alpha = alphas(i2);
        idx = idx + 1;
        paramCell{idx,1} = param;
    end
end
end