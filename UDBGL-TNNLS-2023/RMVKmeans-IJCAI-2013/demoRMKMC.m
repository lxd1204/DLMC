clear;
clc;

load('jaffe_213n_676d_10c_uni.mat')

[resAIO] = RMKMC_grid(Xlow, Xup, Y, nRepeat, mRepeat);
tt = mean(resAIO);
tt(1:4)
%%
load('20NG_3970n_1000d_4c_tfidf_uni.mat');
k = 10;
alpha = 3;
[Xlow, Xup] = IntervalXKnn(X, k, alpha);
inXCell = cell(1,2);
inXCell{1,1} = Xlow;
inXCell{1,2} = Xup;

nRepeat = 10;
mRepeat = 1;
gammaCandi = 10.^[0.1,0.3];
[resAIOCell] = RMKMC_grid(inXCell, y, gammaCandi, nRepeat, mRepeat);
resGridMean = reshape(cell2mat(cellfun(@mean, resAIOCell, 'UniformOutput', false)), ...
    size(resAIOCell{1},2), length(resAIOCell))';
