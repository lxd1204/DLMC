function paramCell = UDBGL_build_param(anchors, alphas, betas)
if ~exist('anchors', 'var')
    anchors = 1000;
end

if ~exist('alphas', 'var')
    alphas = 10.^[-3:3];
end

if ~exist('betas', 'var')
    betas = 10.^[-3:3];
end

nParam = length(anchors) * length(alphas) * length(betas);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(anchors)
    for i2 = 1:length(alphas)
        for i3 = 1:length(betas)
            param = [];
            param.anchor = anchors(i1);
            param.alpha = alphas(i2);
            param.beta = betas(i3);
            idx = idx + 1;
            paramCell{idx,1} = param;
        end
    end
end
end