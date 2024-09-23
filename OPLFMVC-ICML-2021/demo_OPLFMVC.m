prefix = '*';
addpath(genpath(prefix));
data_name = ‘*’；
fprintf(' %s  \n', data_name);
% load data
load([prefix,'datasets\',data_name,'_Kmatrix'],'KH','Y');
% preprocess
numclass = length(unique(Y));
Y(Y==0)=numclass;
Y=double(Y);
numker = size(KH,3);
num = size(KH,1);
KH = kcenter(KH);
KH = knorm(KH);
% algorithm
for irand = 1:50
	s=RandStream('mt19937ar','Seed',irand);
	RandStream.setGlobalStream(s);
	tic;
	[Yout,C,WP,Sigma,obj] = onePassLateFusionMVCBeta(KH,numclass);
	[res_mean(:,irand),res_std(:,irand)]= myNMIACCV2(Yout,Y,numclass);
	timecost(irand) = toc; 
end
% save result
save(savename,'res_mean','res_std','Sigma','timecost' );
clear res_mean res_std KH Y 

