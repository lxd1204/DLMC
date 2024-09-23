function label = TBGL(Xs, nCluster, nAnchor)

nView = length(Xs);
for iView = 1:nView
    a = max(max(Xs{iView}));
    Xs{iView} = double(Xs{iView}./a);
end
%******************************************************
% Parameter Initialization
%******************************************************
alpha = 0.025;                    % the hyper-parameter of L1-norm constraint on E(v)
gama = 0.0001;                    % the hyper-parameter of L12-norm constraint on C(v)
p = 0.9;                          % the p-value of tensor Schatten p-norm
beta = 10;                        % the initialize parameter of Laplacian rank constraint
weight_vector = ones(1, nView)';      % the defult weight_vector of tensor Schatten p-norm

%******************************************************
% Training TBGL
%******************************************************

[C, kesai] = Train_TBGL_v2(Xs, nCluster, nAnchor, weight_vector, p, alpha, beta, gama);
%******************************************************
% Obtain Common Graph
%******************************************************

Common_C = (1/kesai(1))*C{1};
for iView = 2: nView
    Common_C = Common_C + (1/kesai(iView))*C{iView};
end
Sum_kesai = sum(1./kesai);
Common_C = Common_C./Sum_kesai;

%******************************************************
% Obtain clustering results based on  the connectivity of Common_C
%**************************************************
label = my_graphconncomp(Common_C, nCluster, 50);
end