%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep,"tsne_data",filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "FastMICE-TKDE-2023", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_n = 'FMICE';

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
    % FMICE
    %*********************************************************************
    nRepeat = 10;
    ecvis_fmice = zeros(nRepeat, 7);
    time_fmice = zeros(nRepeat, 1);
    fname2 = fullfile(dir_name, [data_name, '_', exp_n, '_res.mat']);
    if ~exist(fname2, 'file')
        for iRepeat = 1:nRepeat
            disp(['FMICE        ', data_name, '        ', num2str(iRepeat), '/', num2str(nRepeat)]);
            fname3 = fullfile(dir_name, [data_name, '_', exp_n, '_repeat', num2str(iRepeat), '.mat']);
            if exist(fname3, 'file')
                clear ecvi t
                load(fname3, 'ecvi', 't');
            else
                t1 = tic;
                [H,label] = runFastMICE(Xs, nCluster);
                t = toc(t1);
                ecvi = ComputeECVIs(Y, label);
                save(fname3, 'ecvi', 't','H','label');
            end
            time_fmice(iRepeat) = t;            
            ecvis_fmice(iRepeat, :) = ecvi';
        end
        FMICE_result_summary = [mean(ecvis_fmice, 1),  mean(time_fmice)];
        save(fname2, 'FMICE_result_summary', 'ecvis_fmice', 'time_fmice');
        disp([data_name, ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear; clc;