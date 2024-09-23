%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep,"tsne_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "LMVSC-AAAI-2020", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'LMVSC';
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
        
        anchor_sizes = [nCluster, 50, 100];
        anchor_sizes = anchor_sizes(anchor_sizes >= nCluster);
        alphas = [0.001, 0.01, 0.1, 1, 10];
        
        paramCell = LMVSC_build_param(anchor_sizes, alphas);
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ecvis_lmvsc_cell = cell(nParam, 1);
        time_lmvsc = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['LMVSC        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear nAnchor alpha bpresFile anchorFile
            
            nAnchor = param.nAnchor;
            alpha = param.alpha;
            
            bpresFile = [file_prefix, '_m', num2str(nAnchor), '_alpha', num2str(alpha), '_res.mat'];
            anchorFile = [file_prefix, '_m', num2str(nAnchor), '.mat'];
            
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
                        [~, Hs{iView}] = litekmeans(Xs{iView}, nAnchor,'MaxIter', 100, 'Replicates', nRepeat_lkm);
                    end
                    t1 = toc(t1_a);
                    save(anchorFile, 'Hs', 't1');
                end
                t2_a = tic;
                [F, label] = lmv(Xs, Y, Hs, alpha);
                t2 = toc(t2_a);
                
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't2', 't1','F','label');
            end
            ecvis_lmvsc_cell{iParam} = ecvi';
            time_lmvsc(iParam) = t1 + t2;
            LMVSC_result_summary = [max(cell2mat(ecvis_lmvsc_cell), [], 1), mean(time_lmvsc)];
            save(dataresFile, 'ecvis_lmvsc_cell', 'time_lmvsc', 'LMVSC_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;