function paramCell = MSGL_build_param(anchor_sizes, alphas, betas, gammas)

if ~exist('anchor_sizes', 'var')
    anchor_sizes = [1:7];
end

if ~exist('alphas', 'var')
    alphas = [1.6];
end


if ~exist('betas', 'var')
    betas = [1.6];
end


if ~exist('gammas', 'var')
    gammas = [1.6];
end

nParam = length(anchor_sizes) * length(alphas) * length(betas) * length(gammas);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(anchor_sizes)
    for i2 = 1:length(alphas)
        for i3 = 1:length(betas)
            for i4 = 1:length(gammas)
                param = [];
                param.anchors = anchor_sizes(i1);
                param.alpha = alphas(i2);
                param.beta = betas(i3);
                param.gamma = gammas(i4);
                idx = idx + 1;
                paramCell{idx,1} = param;
            end
        end
    end
end
end