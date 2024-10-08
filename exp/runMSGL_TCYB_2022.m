%
%
%
clear;
clc;data_path = fullfile(pwd, '..',  filesep,"data_mv", filesep,"tsne_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep,"MBPC_compare_method ",filesep,"MSGL-TCYB-2022", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'MSGL';
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
        
        for iView = 1:nView
            Xs{iView} = Xs{iView}';
        end
        
%         anchor_sizes = [100];
%         anchor_sizes (anchor_sizes < nCluster) = nCluster;
%         anchor_sizes (anchor_sizes > nSmp) = nSmp;
        anchor_sizes = [nCluster];
        alphas = [50];
        % betas = 10.^[-1:1];
        betas = [5*1e4, 1e5, 5*1e5];
        gammas = [-5, -4, -3, -2, -1];
        
        paramCell = MSGL_build_param(alphas, betas, gammas);
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ecvis_msgl_cell = cell(nParam, 1);
        time_msgl = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['MSGL        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear anchors alpha beta gamma bpresFile anchorFile
            
            anchors = param.anchors;
            alpha = param.alpha;
            beta = param.beta;
            gamma = param.gamma;
            bpresFile = [file_prefix, '_m', num2str(anchors), '_alpha', num2str(alpha), '_beta', num2str(beta), '_gamma', num2str(gamma), '_res.mat'];
            anchorFile = [file_prefix, '_m', num2str(anchors), '.mat'];
            
            if exist(bpresFile, 'file')
                clear ecvi t1 t2
                load(bpresFile, 'ecvi', 't1', 't2');
            else
                if exist(anchorFile, 'file')
                    clear Hs t1
                    load(anchorFile, 'Hs', 't1');
                else
                    t1_a = tic;
                    Hs = cell(1, nView);
                    for iView = 1:nView
                        rand('twister',5489);
                        [~, Hi] = litekmeans(Xs{iView}', anchors, 'MaxIter', 100, 'Replicates', nRepeat_lkm);
                        Hs{iView} = Hi';
                    end
                    t1 = toc(t1_a);
                    save(anchorFile, 'Hs', 't1');
                end
                t2_a = tic;
                [label, result,U] = unifiedclusternew_2(Xs, Hs, Y, alpha, beta, gamma, nView);
                t2 = toc(t2_a);
                
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't2', 't1','label','U');
            end
            ecvis_msgl_cell{iParam} = ecvi';
            time_msgl(iParam) = t1 + t2;
            MSGL_result_summary = [max(cell2mat(ecvis_msgl_cell), [], 1), mean(time_msgl)];
            save(dataresFile, 'ecvis_msgl_cell', 'time_msgl', 'MSGL_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;