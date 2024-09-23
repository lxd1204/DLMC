function ComputeBP_mv(Xs, paramCell, file_prefix)

nView = length(Xs);
[nSmp, nFea] = size(Xs{1});

[~, last_folder] = fileparts(file_prefix);

nParam = length(paramCell);
for iParam = 1:nParam
    disp(['ComputeBP        ', last_folder, '        ', num2str(iParam), '/', num2str(nParam)]);
    clear anchor_meta_type anchor_type nAnchor weight_type  beta nNeighbor diffusion_param
    param = paramCell{iParam};
    anchor_meta_type = param.anchor_meta_type;
    nAnchor = param.nAnchor;
    nNeighbor = param.nNeighbor;
    diffusion_param = param.diffusion_param;
    
    if isfield(param, 'anchor_type')
        anchor_type = param.anchor_type;
    end
    
    if isfield(param, 'weight_type')
        weight_type = param.weight_type;
    end
    
    if isfield(param, 'beta')
        beta = param.beta;
    end
    
    
    if exist('anchor_type', 'var') && exist('weight_type', 'var')
        bpdFile =[file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
        bpFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '_', weight_type, '_k', num2str(nNeighbor), '.mat'];
        anchorFile = [file_prefix, '_', anchor_meta_type, '_', anchor_type, '_m', num2str(nAnchor), '.mat'];
    elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
        bpdFile =[file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '_d', num2str(diffusion_param), '.mat'];
        bpFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '_beta', num2str(beta), '_k', num2str(nNeighbor), '.mat'];
        anchorFile = [file_prefix, '_', anchor_meta_type, '_m', num2str(nAnchor), '.mat'];
    end
    
    if ~exist(bpdFile, 'file')
        if exist(bpFile, 'file')
            clear Bs t1 t2 t12;
            if exist('anchor_type', 'var') && exist('weight_type', 'var')
                load(bpFile, 'Bs', 't2', 't1');
            elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                load(bpFile, 'Bs', 't12');
            end
        else
            if exist('anchor_type', 'var') && exist('weight_type', 'var')
                %**************************************
                % Step 1: kmeans for anchor
                %**************************************
                if exist(anchorFile, 'file')
                    clear Xas;
                    load(anchorFile, 'Xas', 't1');
                else
                    tic;
                    if strcmpi(anchor_meta_type, 'hier')
                        nLayer = ceil(log2(nAnchor));
                        Xas = cell(1, nView);
                        for iView = 1:nView
                            X = Xs{iView};
                            switch anchor_type
                                case 'bkhk'
                                    [~, Xa] = hKM(X', (1: nSmp), nLayer, 1);
                                    Xa = Xa';
                                case 'lkm1'
                                    [~, Xa] = litekmeans(X, 2^nLayer, 'Replicates', 1);
                                case 'lkm10'
                                    [~, Xa] = litekmeans(X, 2^nLayer, 'Replicates', 10);
                                case 'kmp'
                                    [~, Xa] = kmeans(X, 2^nLayer, 'Start', 'plus');
                                case 'km'
                                    [~, Xa] = kmeans(X, 2^nLayer, 'Start', 'sample');
                                case 'random'
                                    Xa = X(randperm(nSmp, 2^nLayer), :);
                                otherwise
                                    disp('Not supported yet');
                            end
                            Xas{iView} = Xa;
                        end
                    elseif strcmpi(anchor_meta_type, 'plain')
                        Xas = cell(1, nView);
                        for iView = 1:nView
                            X = Xs{iView};
                            switch anchor_type
                                case 'lkm1'
                                    [~, Xa] = litekmeans(X, nAnchor, 'Replicates', 1);
                                case 'lkm10'
                                    [~, Xa] = litekmeans(X, nAnchor, 'Replicates', 10);
                                case 'kmp'
                                    [~, Xa] = kmeans(X, nAnchor, 'Start', 'plus');
                                case 'km'
                                    [~, Xa] = kmeans(X, nAnchor, 'Start', 'sample');
                                case 'random'
                                    Xa = X(randperm(nSmp, nAnchor), :);
                                otherwise
                                    disp('Not supported yet');
                            end
                            Xas{iView} = Xa;
                        end
                    end
                    t1 = toc;
                    save(anchorFile, 'Xas', 't1');
                end
                
                %**************************************
                % Step 2: construct BP from anchors
                %**************************************
                tic;
                Bs = cell(1, nView);
                for iView = 1:nView
                    X = Xs{iView};
                    Xa = Xas{iView};
                    switch weight_type
                        case 'pkn'
                            B = ConstructBP_pkn(X, Xa, 'nNeighbor', nNeighbor);
                        case 'lkr'
                            B = ConstructBP_lkr(X, Xa, 'nNeighbor', nNeighbor);
                    end
                    Bs{iView} = B;
                end
                
                t2 = toc;
                save(bpFile, 'Bs', 't2', 't1');
            elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                tic;
                Bs = cell(1, nView);
                for iView = 1:nView
                    X = Xs{iView};
                    B = ConstructionBP_onestep(X, nAnchor, beta);
                    Bs{iView} = B;
                end
                t12 = toc;
                save(bpFile, 'Bs', 't12');
            end
        end
        %**************************************
        % Step 3: diffusion and sparse BP
        %**************************************
        if diffusion_param > 0
            tic;
            for iView = 1:nView
                Bs{iView} = DiffusionBP_ppr(Bs{iView}, 'alpha', diffusion_param, 'nNeighbor', nNeighbor);
            end
            t3 = toc;
            if exist('anchor_type', 'var') && exist('weight_type', 'var')
                save(bpdFile, 'Bs', 't3', 't2', 't1');
            elseif exist('beta', 'var') && ~exist('anchor_type', 'var') && ~exist('weight_type', 'var')
                save(bpdFile, 'Bs', 't3', 't12');
            end
        end
    end
end

end