%
%
%
clear;
clc;
% data_path = fullfile(pwd, '..',  filesep, "data_mv",filesep,"data_mv_all_format",filesep);
data_path = fullfile(pwd, '..',  filesep, "data_mv",filesep,"MV-LargeScale-Data",filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "SBP", filesep);
addpath(code_path);
code_path = fullfile(pwd, '..',  filesep, "Ours", filesep);
addpath(code_path);

dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

% exp_name = 'EXP_result_try_square1_knn20';
exp_name = 'EXP_result_try_0.5';
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
        disp(data_name);
        load(data_name);
        Xs = data;
        Y = truth;
        nView = length(Xs);
        for iView = 1:nView
            Xs{iView} = double(full(Xs{iView}));
        end
        nSmp = size(Xs{1}, 1);
        nCluster = length(unique(Y));
        if nSmp < 1000000
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
                % weight_types = {'pkn', 'lkr'};
                weight_types = {'pkn'};
                % knn_sizes = [5, 10, 15, 20];
                knn_sizes = [5];
                % diffusion_params = [0, 3, 5, 7, 9];
                diffusion_params = [5];
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
                anchor_types = {'lkm1'};
                weight_types = {'pkn'};
%                 knn_sizes = [5, 10, 15, 20,30,40,50];
                knn_sizes = [20];
                % diffusion_params = [0, 3, 5, 7, 9];
                % diffusion_params = [0.05, 0.1, 0.15, 0.9];
                diffusion_params = [0];
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
                knn_sizes = [5, 10, 15, 20];
                % knn_sizes = [5, 10];
                diffusion_params = [0, 3, 5, 7, 9];
                % diffusion_params = [0, 3];
                betas = [0.1, 1, 10];
                paramCell = ConstructBP_build_param(anchor_meta_type, anchor_types, anchor_sizes, weight_types, knn_sizes, diffusion_params, betas);
            end
            nParam = length(paramCell);
            
            nRepeat_lkm = 10;
            
            ComputeBP_mv(Xs, paramCell, file_prefix);
            
            ecvis_ourv4mvablation_cell = cell(nParam, 1);
            time_ourv4mvablation = zeros(nParam, 1);
            param_idx = [];
            p2s_old = [];
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
                    bpresFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_res.mat'];
                    bpFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '.mat'];
                    anchorFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '.mat'];
                elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                    bpresFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '_res.mat'];
                    bpFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '.mat'];
                    anchorFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '.mat'];
                end
                
                
                if exist(bpresFile, 'file')
                    clear ecvi t4 tall
                    load(bpresFile, 'ecvi', 't_all', 't4');
                else
                    % bpdFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
                    clear Bs t1 t2 t12;
                    if exist('anchor_type', 'var') && exist('weight_type', 'var')
                        load(bpFile, 'Bs', 't2', 't1');
                    elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                        load(bpFile, 'Bs', 't12');
                    end
                    tic;
                    [label, U, Vs, objHistory,Y_normalized] = MBPLR_try(Bs, nCluster);
                    plot(objHistory);
                    t4 = toc;
                    if exist('anchor_type', 'var') && exist('weight_type', 'var')
                        t_all = t4 + t1 + t2;
                    elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                        t_all = t4 + t12;
                    end
                    ecvi = ComputeECVIs(Y, label)';
                    save(bpresFile, 'label', 'ecvi', 't_all', 't4','Y_normalized');
                end                
                ecvis_ourv4mvablation_cell{iParam} = ecvi(:)';
                time_ourv4mvablation(iParam) = t_all;
                
                p2s_new = [anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor)];
                if ~strcmpi(p2s_old, p2s_new)
                    param_idx = [param_idx; iParam]; %#ok
                    p2s_old = p2s_new;
                end
            end
            ecvis_ourv4mvablation_table = cell2mat(ecvis_ourv4mvablation_cell);
            ecvis_ourv4mvablation_table = ecvis_ourv4mvablation_table(param_idx, :);
            time_ourv4mvablation = time_ourv4mvablation(param_idx);
            Ourv4mvablation_result_summary = [max(ecvis_ourv4mvablation_table, [], 1), mean(time_ourv4mvablation)];
            save(dataresFile, 'ecvis_ourv4mvablation_cell', 'ecvis_ourv4mvablation_table', 'paramCell', 'Ourv4mvablation_result_summary' ,'time_ourv4mvablation');
        end
    end 
end % end of data sets


rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;