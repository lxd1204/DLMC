%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "USENC_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "data_mv_standard", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "USENC-TKDE-2020", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_n = 'USENC';

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
    % USENC
    %*********************************************************************
    nRepeat = 10;
    m = 20; % Ensemble size
    distance = 'euclidean';
    p = 1000;
    Knn = 5;
    bcsLowK = nCluster;
    bcsUpK = 2 * nCluster;
    ecvis_usenc = zeros(nRepeat, 7);
    time_usenc = zeros(nRepeat, 1);
    fname2 = fullfile(dir_name, [data_name, '_', exp_n, '_res.mat']);
    if ~exist(fname2, 'file')
        for iRepeat = 1:nRepeat
            disp(['USENC        ', data_name, '        ', num2str(iRepeat), '/', num2str(nRepeat)]);
            fname3 = fullfile(dir_name, [data_name, '_', exp_n, '_repeat', num2str(iRepeat), '.mat']);
            if exist(fname3, 'file')
                clear ecvi t
                load(fname3, 'ecvi', 't');
            else
                t1 = tic;
                [H,label] = USENC_mv(Xs, nCluster);
                t = toc(t1);
                ecvi = ComputeECVIs(Y, label);
                save(fname3, 'ecvi', 't','H','label');
             end
            time_usenc(iRepeat) = t;
            ecvis_usenc(iRepeat, :) = ecvi';
        end
        USENC_result_summary = [mean(ecvis_usenc, 1),  mean(time_usenc)];
        save(fname2, 'USENC_result_summary', 'ecvis_usenc', 'time_usenc');
        disp([data_name, ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear; clc;