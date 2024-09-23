%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "tsne_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "SMCMB-TBD-2023", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'SMCMB';
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
            Xs{iView} = double(full(Xs{iView}));
        end
        
        anchor_sizes = [nCluster, 50, 100];
        anchor_sizes (anchor_sizes < nCluster) = nCluster;
        anchor_sizes (anchor_sizes > nSmp) = nSmp;
        
        ms = 5;
        betas = 10.^[-5:5];
        paramCell = SMCMB_build_param(anchor_sizes, ms, betas);
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ecvis_smcmb_cell = cell(nParam, 1);
        time_smcmb = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['SMCMB        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear nAnchor  bpresFile
            
            nAnchor = param.anchor;
            m = param.m;
            beta = param.beta;
            
            bpresFile = [file_prefix, '_nAnchor', num2str(nAnchor), '_m', num2str(m), '_beta', num2str(beta), '_res.mat'];
            
            if exist(bpresFile, 'file')
                clear ecvi
                load(bpresFile, 'ecvi', 't');
            else
                t1_a = tic;
                [H,label] = SMCMB(Xs, nCluster, m, nAnchor, beta);
                t = toc(t1_a);
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't','H','label');
            end
            ecvis_smcmb_cell{iParam} = ecvi';
            time_smcmb(iParam) = t;
            SMCMB_result_summary = [max(cell2mat(ecvis_smcmb_cell), [], 1), mean(time_smcmb)];
            save(dataresFile, 'ecvis_smcmb_cell', 'time_smcmb', 'SMCMB_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
% clear;clc;