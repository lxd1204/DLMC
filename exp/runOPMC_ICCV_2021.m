%
%
%
clear;
clc;
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "data_mv_standard", filesep);
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "tsne_data", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "OPMC-ICCV-2021", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_n = 'OPMC';

for i1 = 1:length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = [pwd, filesep, exp_n, filesep, data_name];
    create_dir(dir_name);
    
    clear Xs Y;
    load(data_name);
    nView = length(Xs);
    nSmp = size(Xs{1}, 1);
    nCluster = length(unique(Y));
    
    %*********************************************************************
    % OPMC
    %*********************************************************************
    nRepeat = 10;
    ecvis_opmc = zeros(nRepeat, 7);
    time_opmc = zeros(nRepeat, 1);
    error_opmc = zeros(nRepeat, 1);
    fname2 = fullfile(dir_name, [data_name, '_', exp_n, '_res.mat']);
    if ~exist(fname2, 'file')
        for iRepeat = 1:nRepeat
            
            disp(['OPMC        ', data_name, '        ', num2str(iRepeat), '/', num2str(nRepeat)]);
            fname3 = fullfile(dir_name, [data_name, '_', exp_n, '_repeat', num2str(iRepeat), '.mat']);
            if exist(fname3, 'file')
                clear ecvi t error
                load(fname3, 'ecvi', 't', 'error');
            else
                t1 = tic;
                [label, C, W, beta, obj] = opmc(Xs, nCluster);
                t = toc(t1);
                error = obj(end);
                ecvi = ComputeECVIs(Y, label);
                save(fname3, 'ecvi', 't', 'error');
            end
            time_opmc(iRepeat) = t;
            error_opmc(iRepeat) = error;
            ecvis_opmc(iRepeat, :) = ecvi';
        end
        OPMC_result_summary = [mean(ecvis_opmc, 1),  mean(time_opmc)];
        save(fname2, 'OPMC_result_summary', 'ecvis_opmc', 'time_opmc', 'error_opmc');
        disp([data_name, ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear; clc;