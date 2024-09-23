%
%
%
clear;
clc;
data_path = fullfile(pwd,  '..',  filesep, "data_mv", filesep,"tsne_data",filesep);
% data_path = fullfile(pwd, '..',  filesep, "data_mv", filesep, "clean", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = fullfile(pwd, '..',  filesep, "FDAGF-AAAI-2023", filesep);
addpath(genpath(code_path));


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};

exp_name = 'FDAGF';
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
        
        anchor_traverse = 4; % R
        anchor_sizes = (1:anchor_traverse)*nCluster;
        anchor_sizes (anchor_sizes < nCluster) = nCluster;
        anchor_sizes (anchor_sizes > nSmp) = nSmp;

        alphas = 10.^[-5, -1, 1, 3];
        lambdas = 10.^[1, 3, 5];
        
        paramCell = FDAGF_build_param(alphas, lambdas);
        nParam = length(paramCell);
        
        nRepeat_lkm = 1;
        
        ecvis_fdagf_cell = cell(nParam, 1);
        time_fdagf = zeros(nParam, 1);
        for iParam = 1:nParam
            disp(['FDAGF        ', data_name, '        ', num2str(iParam), '/', num2str(nParam)]);
            param = paramCell{iParam};
            clear anchors alpha beta lambda bpresFile anchorFile
            
            anchors = anchor_sizes;
            alpha = param.alpha;
            lambda = param.lambda;
            str_a = arrayfun(@num2str, anchors, 'UniformOutput', false);
            anchor_str = strjoin(str_a, '_');
            
            bpresFile = [file_prefix, '_m', anchor_str, '_alpha', num2str(alpha), '_lambda', num2str(lambda), '_res.mat'];
            anchorFile = [file_prefix, '_m', anchor_str, '.mat'];
            
            if exist(bpresFile, 'file')
                clear ecvi t1 t2
                load(bpresFile, 'ecvi', 't1', 't2');
            else
                if exist(anchorFile, 'file')
                    clear Hs t1
                    load(anchorFile, 'Hs', 't1');
                else
                    t1_a = tic;
                    Hs = cell(nView, length(anchors));
                    for iView = 1:nView
                        for iAnchor = 1:length(anchors)
                            rand('twister',5489);
                            [~, Hi] = litekmeans(Xs{iView}', anchors(iAnchor), 'MaxIter', 100, 'Replicates', nRepeat_lkm);
                            Hs{iView, iAnchor} = Hi';
                        end
                    end
                    t1 = toc(t1_a);
                    save(anchorFile, 'Hs', 't1');
                end
                t2_a = tic;
                [label,U] = algorithm(Xs, Y, Hs, alpha, lambda, anchors);
                t2 = toc(t2_a);
                
                ecvi = ComputeECVIs(Y, label);
                save(bpresFile, 'ecvi', 't2', 't1','label','U');
            end
            ecvis_fdagf_cell{iParam} = ecvi';
            time_fdagf(iParam) = t1 + t2;
            FDAGF_result_summary = [max(cell2mat(ecvis_fdagf_cell), [], 1), mean(time_fdagf)];
            save(dataresFile, 'ecvis_fdagf_cell', 'time_fdagf', 'FDAGF_result_summary');
        end
        disp([data_name(1:end-4), ' has been completed!']);
    end
end

rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);
% clear;clc;