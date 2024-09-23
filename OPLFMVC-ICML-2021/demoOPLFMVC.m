
dataName ={'CSTR_476n_1000d_4c_uni_12k',...
    'Prostate_GE_102n_5966d_2c_12k'...
    'RELATHE_1427n_4322d_2c_tfidf_uni_12k',...
    'USPS49_1673n_256d_2c_12k',...
    'ucidigit_2000n_6v_12k_10c',... 
    };
for i1 =1:length(dataName)
    load(fullfile(pwd, '..', 'data', dataName{i1}));
    nRepeat = 10;
    alphas = 10.^[-4:1:1]; % trade-off
    betas = [1,5,10,15,20,25,30]; % local regularization
    gammas = [0.1,1,5,10,15,20,25,30];
    Ks2 = zeros(size(Ks));
    for i2 = 1:size(Ks,3)
        Ktmp = knorm(Ks(:,:,i2));
        if max(Ktmp(:)) > 1
            Ktmp = Ktmp  / max(Ktmp(:));
        end
        Ks2(:,:,i2) = Ktmp;
    end
    clear Ks;
    [bstResult1, timeElapsed1, res_oplfmvc_aio] = OPLFMVC_single_dataset(Ks2, Y, nRepeat);
    save([dataName{i1}, '_oplfmvc_multi_kernel.mat'],  'res_oplfmvc_aio','timeElapsed1','bstResult1');
end