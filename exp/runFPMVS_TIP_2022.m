%
%
%
clear;
clc;
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep,"tsne_data", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean_v2", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "FPMVS-TIP-2022", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_n = 'FPMVS';

for i1 = 1:length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = [pwd, filesep, exp_n, filesep, data_name];
    create_dir(dir_name);
    
    clear Xs Y;
    load(data_name);
    nView = length(Xs);
    for iView = 1:nView
        Xs{iView} = double(Xs{iView});
    end
    nSmp = size(Xs{1}, 1);
    nCluster = length(unique(Y));
    if nSmp < 20000
        %*********************************************************************
        % FPMVS
        %*********************************************************************
        nRepeat = 10;
        %*********************************************************************
        % Default setting as the paper, check the result of Caltech101-20
        %*********************************************************************
        anchor = nCluster ;
        d = (1)*nCluster ;
        lambda=0;
        ecvis_fpmvs = zeros(nRepeat, 6);
        time_fpmvs = zeros(nRepeat, 1);
        fname2 = fullfile(dir_name, [data_name, '_', exp_n, '_res.mat']);
        if ~exist(fname2, 'file')
            disp(['FPMVS        ', data_name, '        Consensus Graph Learning']);
            t1 = tic;
            [U, A, W, Z, iter, obj, alpha, P] = algo_qp(Xs, Y, lambda, d, anchor); % X,Y,lambda,d,numanchor
            t2 = toc(t1);
            stream = RandStream.getGlobalStream;
            reset(stream);
            U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,size(U,2));
            maxIter = 50;
            
            for iRepeat = 1:nRepeat
                disp(['FPMVS        ', data_name, '        ', num2str(iRepeat), '/', num2str(nRepeat)]);
                fname3 = fullfile(dir_name, [data_name, '_', exp_n, '_repeat', num2str(iRepeat), '.mat']);
                if exist(fname3, 'file')
                    clear ecvi t4
                    load(fname3, 'ecvi', 't4');
                else
                    t3 = tic;
                    label = litekmeans(U_normalized, nCluster, 'MaxIter',100, 'Replicates',1);
                    ecvi = ComputeECVIs(Y, label);
                    t4 = toc(t3);
                    save(fname3, 'ecvi', 't4','U_normalized','label');
                end
                ecvis_fpmvs(iRepeat, :) = ecvi';
                time_fpmvs(iRepeat) = t2 + t4;
            end
            
            FPMVS_result_summary = [mean(ecvis_fpmvs, 1),  mean(time_fpmvs)];
            save(fname2, 'FPMVS_result_summary', 'ecvis_fpmvs', 'time_fpmvs');
            disp([data_name, ' has been completed!']);
        end
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear; clc;