function paramCell = BPParamConfiguration_plain(anchor_types, anchor_sizes, weight_types, diffusion_params)

%**********************************************
% Step 1: prepare grid of BP
%**********************************************
if ~exist('anchor_types', 'var')
    anchor_types = {'lkm1', 'lkm10', 'km', 'kmp', 'random'};
end
if ~exist('anchor_sizes', 'var')
    anchor_sizes = [5];
end

if ~exist('weight_types', 'var')
    weight_types = {'pkn', 'lkr'};
end

if ~exist('knn_sizes', 'var')
    knn_sizes = [5, 10, 15, 20];
end

if ~exist('diffusion_params', 'var')
    diffusion_params = [0, 3, 5, 7, 9];
end

anchor_meta_type = 'plain';

paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params);
end