function paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params, betas)

if ~exist('anchor_types', 'var')
    anchor_types = {'km', 'km+', 'random'};
end

if ~exist('anchor_sizes', 'var')
    anchor_sizes = 100;
end

if ~exist('weight_types', 'var')
    weight_types = {'pkn', 'lkr'};
end

if ~exist('betas', 'var')
    betas = [0.1, 1, 10];
end

if ~exist('knn_sizes', 'var')
    knn_sizes = [5, 10, 15, 20];
end

if ~exist('diffusion_params', 'var')
    diffusion_params = [0, 3, 5, 7, 9];
end

if strcmpi(anchor_meta_type, 'hier') || strcmpi(anchor_meta_type, 'plain')
    nParam = length(anchor_types) * length(anchor_sizes) * length(weight_types) * length(knn_sizes) * length(diffusion_params);
    paramCell = cell(nParam, 1);
    idx = 0;
    for i1 = 1:length(anchor_types)
        for i2 = 1:length(anchor_sizes)
            for i3 = 1:length(weight_types)
                for i4 = 1:length(knn_sizes)
                    for i5 = 1:length(diffusion_params)
                        param = [];
                        param.anchor_meta_type = anchor_meta_type;
                        param.anchor_type = anchor_types{i1};
                        param.nAnchor = anchor_sizes(i2);
                        param.weight_type = weight_types{i3};
                        param.nNeighbor = knn_sizes(i4);
                        param.diffusion_param = diffusion_params(i5);
                        idx = idx + 1;
                        paramCell{idx,1} = param;
                    end
                end
            end
        end
    end
elseif strcmpi(anchor_meta_type, 'onestep')
    nParam = length(anchor_sizes) * length(betas) * length(knn_sizes) * length(diffusion_params);
    paramCell = cell(nParam, 1);
    idx = 0;
    for i1 = 1:length(anchor_sizes)
        for i2 = 1:length(betas)
            for i3 = 1:length(knn_sizes)
                for i4 = 1:length(diffusion_params)
                    param = [];
                    param.anchor_meta_type = anchor_meta_type;
                    param.nAnchor = anchor_sizes(i1);
                    param.beta = betas(i2);
                    param.nNeighbor = knn_sizes(i3);
                    param.diffusion_param = diffusion_params(i4);
                    idx = idx + 1;
                    paramCell{idx,1} = param;
                end
            end
        end
    end
end
end