%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep,  "data_mv", filesep, "tsne_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "EMKMC-TKDE-2023", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'EMKMC';
for i1 = 1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    base_dir_name = [pwd, filesep,  exp_name, filesep, data_name];
    create_dir(base_dir_name);
    file_prefix = [base_dir_name, filesep, data_name];
    data_full_path = strcat(data_path, datasetCandi{i1});
    
    dataresFile = [file_prefix, '_', exp_name, '_res.mat'];
    if ~exist(dataresFile, 'file')
        clear Xs Y;
        load(data_name);
        nView = length(Xs);
        nSmp = size(Xs{1}, 1);
        nCluster = length(unique(Y));
        
        
        anchor_sizes = [nCluster * (1 : 6)'; 1000] * ones(1, nView);  % m * nViews
        anchor_sizes (anchor_sizes < nCluster) = nCluster;
        anchor_sizes (anchor_sizes > nSmp) = nSmp;
        % anchor_sizes =[10, 16, 14, 12, 12, 16];
        gammas = [1.5];
        
        paramCell = EMKMC_build_param(anchor_sizes, gammas);
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ecvis_emkmc_cell = cell(nParam, 1);
        time_emkmc = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['EMKMC        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear anchors gamma bpresFile anchorFile
            
            anchors = param.anchors;
            gamma = param.gamma;
            str_a = arrayfun(@num2str, anchors, 'UniformOutput', false);
            anchor_str = strjoin(str_a, '_');
            
            bpresFile = [file_prefix, '_m', anchor_str, '_gamma', num2str(gamma), '_res.mat'];
            bpFile = [file_prefix, '_m', anchor_str, '.mat'];
            
            if exist(bpresFile, 'file')
                clear ecvi t1 t2
                load(bpresFile, 'ecvi', 't1', 't2');
            else
                if exist(bpFile, 'file')
                    clear Bs t1
                    load(bpFile, 'Bs', 't1');
                else
                    t1_a = tic;
                    Bs = cell(1, nView);
                    for iView = 1:nView
                        rand('twister',5489);
                        [~, Anchors] = litekmeans(Xs{iView}, anchors(iView), 'MaxIter', 100, 'Replicates', nRepeat_lkm);
                        Bi = ConstructA_NP(Xs{iView}', Anchors');
                        Bs{iView} = Bi;
                    end
                    t1 = toc(t1_a);
                    save(bpFile, 'Bs', 't1');
                end
                t2_a = tic;
                [F,label] = EMKMC_Bs(Bs, nCluster, gamma);
                t2 = toc(t2_a);
                
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't2', 't1','F','label');
            end
            ecvis_emkmc_cell{iParam} = ecvi';
            time_emkmc(iParam) = t1 + t2;
            EMKMC_result_summary = [max(cell2mat(ecvis_emkmc_cell), [], 1), mean(time_emkmc)];
            save(dataresFile, 'ecvis_emkmc_cell', 'time_emkmc', 'EMKMC_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;