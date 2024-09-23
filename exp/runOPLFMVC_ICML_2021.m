%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "tsne_data", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "SBP", filesep);
addpath(code_path);
code_path = fullfile(pwd, '..',  filesep, "OPLFMVC-ICML-2021", filesep);
addpath(code_path);

dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'OPLFMVC';
res_aio = cell(length(datasetCandi), 4);
time_aio = cell(length(datasetCandi), 4);
for i1 = 1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    base_dir_name = [pwd, filesep,  exp_name, filesep, data_name];
    create_dir(base_dir_name);
    file_prefix = [base_dir_name, filesep, data_name];
    data_full_path = strcat(data_path, datasetCandi{i1});
    
    anchor_meta_type = 'plain';
    dataresFile = [file_prefix, '_', anchor_meta_type, '_', exp_name, '_res.mat'];
    if ~exist(dataresFile, 'file')
        clear Xs Y;
        load(data_name);
        nView = length(Xs);
        nSmp = size(Xs{1}, 1);
        nCluster = length(unique(Y));
        
        
        if strcmpi(anchor_meta_type, 'hier')
            if nSmp > 50000
                anchor_sizes = [512, 1024, 2048, 4096];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 10000
                anchor_sizes = [256, 512, 1024, 2048];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 5000
                anchor_sizes = [128, 256, 512, 1024];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            else
                anchor_sizes = [64, 128, 256, 512];
                anchor_sizes = anchor_sizes(anchor_sizes < nSmp);
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            end
            % anchor_types = {'bkhk', 'lkm1', 'lkm10', 'km', 'kmp', 'random'};
            anchor_types = {'bkhk'};
            weight_types = {'pkn'};
            % knn_sizes = [5, 10, 15, 20];
            knn_sizes = [5, 10];
            % diffusion_params = [0, 3, 5, 7, 9];
            diffusion_params = [0, 3];
            paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params);
        elseif strcmpi(anchor_meta_type, 'plain')
            if nSmp > 50000
                anchor_sizes = [1000, 1500, 2000, 2500, 3000, 3500, 4000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 10000
                anchor_sizes = [500, 800, 1000, 1200, 1500, 2000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 5000
                anchor_sizes = [100, 200, 500, 800,1000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            else
                anchor_sizes = [100, 200, 300, 400, 500];
                anchor_sizes = anchor_sizes(anchor_sizes < nSmp);
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            end
            % anchor_types = {'lkm1', 'lkm10'};
            anchor_types = {'lkm10'};
            weight_types = {'pkn'};
            knn_sizes = [5, 10];
            % knn_sizes = [5, 10];
            diffusion_params = [0, 3];
            % diffusion_params = [0, 3];            
            paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params);
        elseif strcmpi(anchor_meta_type, 'onestep')
            if nSmp > 50000
                anchor_sizes = [1000, 1500, 2000, 2500, 3000, 3500, 4000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 10000
                anchor_sizes = [500, 800, 1000, 1200, 1500, 2000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            elseif nSmp > 5000
                anchor_sizes = [100, 200, 500, 800,1000];
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            else
                anchor_sizes = [100, 200, 300, 400, 500];
                anchor_sizes = anchor_sizes(anchor_sizes < nSmp);
                anchor_sizes = anchor_sizes(anchor_sizes > nCluster);
            end
            anchor_types = [];
            weight_types = [];
            % knn_sizes = [5, 10, 15, 20];
            knn_sizes = [5, 10];
            % diffusion_params = [0, 3, 5, 7, 9];
            diffusion_params = [0, 3];
            betas = [0.1, 1, 10];
            paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params, betas);
        end
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ComputeBP_mv(Xs, paramCell, file_prefix);
        
        ecvis_oplfmvc_cell_a = cell(nParam, 1);
        ecvis_oplfmvc_cell_b = cell(nParam, 1);
        time_oplfmvca = zeros(nParam, 1);
        time_oplfmvcb = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['CutBP        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear anchor_meta_type anchor_type nAnchor weight_type  beta nNeighbor diffusion_param
            anchor_meta_type = param.anchor_meta_type;
            nAnchor = param.nAnchor;
            nNeighbor = param.nNeighbor;
            diffusion_param = param.diffusion_param;
            
            if isfield(param, 'anchor_type')
                anchor_type = param.anchor_type;
            end
            
            if isfield(param, 'weight_type')
                weight_type = param.weight_type;
            end
            
            if isfield(param, 'beta')
                beta = param.beta;
            end
            
            
            if exist('anchor_type', 'var') && exist('weight_type', 'var')
                bpresFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '_res.mat'];
                bpdFile =[file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
                bpFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '.mat'];
                anchorFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '.mat'];
            elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                bpresFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '_res.mat'];
                bpdFile =[file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
                bpFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '.mat'];
                anchorFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '.mat'];
            end
            
            
            if exist(bpresFile, 'file')
                clear ecvi_a ecvi_b t4_a t4_b
                load(bpresFile, 'ecvi_a', 'ecvi_b', 't4_a', 't4_b');
            else
                % bpdFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
                clear Bs t1 t2 t12;
                if exist('anchor_type', 'var') && exist('weight_type', 'var')
                    load(bpFile, 'Bs', 't2', 't1');
                elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                    load(bpFile, 'Bs', 't12');
                end
                tic;
                Us_a = zeros(nSmp, nCluster, nView);
                for iView = 1:nView
                    Us_a(:, :, iView) = bcut_lsc(Bs{iView}, nCluster);
                end
                WP = zeros(nCluster, nCluster, nView);
                for iView = 1:nView
                    WP(:, :, iView) = eye(nCluster);
                end                
                Uc = reshape(Us_a, nSmp, nCluster * nView);
                Uc =NormalizeFea(Uc);
                [~, C] = kmeans(Uc, nCluster, 'MaxIter', 100, 'Replicates', nRepeat_lkm);
                C = orth(C);
                [Y_out, ~, ~, ~, ~] = OPLFMVC(Us_a, WP, C, nCluster);
                [~, label_a] = max(Y_out, [], 2);
                t4_a = toc;
                ecvi_a = ComputeECVIs(Y, label_a)';
                
                tic;
                Us_b = zeros(nSmp, nCluster, nView);
                for iView = 1:nView
                    Us_b(:, :, iView) = bcut_tcut(Bs{iView}, nCluster);
                end
                WP = zeros(nCluster, nCluster, nView);
                for iView = 1:nView
                    WP(:, :, iView) = eye(nCluster);
                end                
                Uc = reshape(Us_b, nSmp, nCluster * nView);
                Uc =NormalizeFea(Uc);
                [~, C] = kmeans(Uc, nCluster, 'MaxIter', 100, 'Replicates', nRepeat_lkm);
                C = orth(C);
                [Y_out, ~, ~, ~, ~] = OPLFMVC(Us_b, WP, C, nCluster);
                [~, label_b] = max(Y_out, [], 2);
                t4_b = toc;
                ecvi_b = ComputeECVIs(Y, label_b)';
                if exist('anchor_type', 'var') && exist('weight_type', 'var')
                    t4_a = t4_a + t1 + t2;
                    t4_b = t4_b + t1 + t2;
                elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                    t4_a = t4_a + t12;
                    t4_b = t4_b + t12;
                end
                save(bpresFile, 'ecvi_a', 'ecvi_b', 't4_a', 't4_b','Y_out','label_b');
            end
            ecvis_oplfmvc_cell_a{iParam} = ecvi_a;
            ecvis_oplfmvc_cell_b{iParam} = ecvi_b(:)';
            time_oplfmvca(iParam) = t4_a;
            time_oplfmvcb(iParam) = t4_b;
        end
        ecvis_oplfmvc_table_a = cell2mat(ecvis_oplfmvc_cell_a);
        ecvis_oplfmvc_table_b = cell2mat(ecvis_oplfmvc_cell_b);
        OPLFMVCa_result_summary = [max(ecvis_oplfmvc_table_a, [], 1), mean(time_oplfmvca)];
        OPLFMVCb_result_summary = [max(ecvis_oplfmvc_table_b, [], 1), mean(time_oplfmvcb)];
        save(dataresFile, 'ecvis_oplfmvc_cell_a', 'ecvis_oplfmvc_cell_b', 'ecvis_oplfmvc_table_a', 'ecvis_oplfmvc_table_b', 'paramCell', 'OPLFMVCa_result_summary', 'OPLFMVCb_result_summary', 'time_oplfmvca', 'time_oplfmvcb');
    end
end % end of data sets


rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;