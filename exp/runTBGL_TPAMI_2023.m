%
%
%
clear;
clc;
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "data_emkmc", filesep);
data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "data_tiny", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "TBGL-TPAMI-2023", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'TBGL';
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
        if nSmp < 10000

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
            knn_size = 10;
            anchor_sizes(anchor_sizes  < knn_size + 1) = [];
            paramCell = TBGL_build_param(anchor_sizes);
            nParam = length(paramCell);

            % nRepeat_lkm = 10;

            ecvis_sfmc_cell = cell(nParam, 1);
            time_sfmc = zeros(nParam, 1);
            for iParam = 1:nParam
                disp(['TBGL        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
                param = paramCell{iParam};
                clear nAnchor  bpresFile

                nAnchor = param.nAnchor;

                bpresFile = [file_prefix, '_nAnchor', num2str(nAnchor), '_res.mat'];

                if exist(bpresFile, 'file')
                    clear ecvi
                    load(bpresFile, 'ecvi', 't');
                else
                    t1_a = tic;
                    label = TBGL(Xs, nCluster, nAnchor);
                    t = toc(t1_a);
                    ecvi = ComputeECVIs(Y, label);
                    save(bpresFile, 'ecvi', 't');
                end
                ecvis_sfmc_cell{iParam} = ecvi';
                time_sfmc(iParam) = t;
                TBGL_result_summary = [max(cell2mat(ecvis_sfmc_cell), [], 1), mean(time_sfmc)];
                save(dataresFile, 'ecvis_sfmc_cell', 'time_sfmc', 'TBGL_result_summary');
            end
            disp([data_name(1:end-4), ' has been completed!']);
        end
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
clear;clc;