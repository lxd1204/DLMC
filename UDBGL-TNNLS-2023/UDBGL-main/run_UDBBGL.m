%
%
%
clear;
close all;
clc;
data_path = fullfile(pwd, filesep, "MBPC_data", filesep);
addpath(data_path);

dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};
for i1 = 1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    load(data_name);
    disp(data_name);
    
    result_path = fullfile(pwd,filesep,'MBPC_data_result',filesep,[data_name,'_result']);
    fid = fopen(result_path,'a');
    X = Xs;
    n = numel(Y);
    c = numel(unique(Y)); % The number of clusters
    
    % Parameters
    m = c; %The number of anchors
%     alpha = 1e-3;
%     beta = 1e-5;  
    
    alpha = 10.^(-3:3); %对照论文对参数就进行了改进
    beta = 10.^(-3:3);  %对照论文对参数就进行了改进
    for  i = 1:length(alpha)
        for  j = 1:length(beta)
            tic;
            Label = UDBGL(X,c,m,alpha(i),beta(j));
            time = toc;
            result = my_eval_y(Label, Y);
            fprintf('Dataset:%s\t beta:%.4f\t  alpha:%.4f\t ACC:%.4f\t NMI:%.4f\t Purity:%.4f\t ARI:%.4f\t RI:%.4f\t  MI:%.4f\t HI:%.4f\t Fscore:%.4f\t Precision:%.4f\t Recall:%.4f\t Entropy:%.4f\tSDCS:%.4f\t RME:%.4f\t Time:%.4f\n',data_name,alpha(i),beta(j),result',time);
            fprintf(fid,'Dataset:%s\t beta:%.4f\t  alpha:%.4f\t ACC:%.4f\t NMI:%.4f\t Purity:%.4f\t ARI:%.4f\t RI:%.4f\t  MI:%.4f\t HI:%.4f\t Fscore:%.4f\t Precision:%.4f\t Recall:%.4f\t Entropy:%.4f\tSDCS:%.4f\t RME:%.4f\t Time:%.4f\n',data_name,alpha(i),beta(j), result',time);
        end
    end
 end