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
code_path = fullfile(pwd, '..',  filesep, "EOMSCCA-AAAI-2022", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'EOMSC';
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
        
        anchor_sizes = nCluster * [1 2 3 4 5 6 7]; 
        anchor_sizes (anchor_sizes < nCluster) = nCluster;
        anchor_sizes (anchor_sizes > nSmp) = nSmp;
        d_rates = nCluster * [1 2 3 4 5 6 7];
        
        paramCell = EOMSC_build_param(anchor_sizes, d_rates);
        nParam = length(paramCell);
        
        nRepeat_lkm = 10;
        
        ecvis_eomsc_cell = cell(nParam, 1);
        time_eomsc = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['EOMSC        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear nAnchor d_rate bpresFile anchorFile
            
            anchors = param.anchors;
            d_rate = param.d_rate;
            
            bpresFile = [file_prefix, '_m', num2str(anchors), '_param_d', num2str(d_rate), '_res.mat'];
            bpFile = [file_prefix, '_m', num2str(anchors), '.mat'];
            
            if exist(bpresFile, 'file')
                clear ecvi 
                load(bpresFile, 'ecvi', 't');
            else
                t1_a = tic;
                [A,W,Z,iter,obj,alpha,label] = algo_qp(Xs, Y, d_rate, anchors);
                t = toc(t1_a);
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't');
            end
            ecvis_eomsc_cell{iParam} = ecvi';
            time_eomsc(iParam) = t;
            EOMSC_result_summary = [max(cell2mat(ecvis_eomsc_cell), [], 1), mean(time_eomsc)];
            save(dataresFile, 'ecvis_eomsc_cell', 'time_eomsc', 'EOMSC_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;