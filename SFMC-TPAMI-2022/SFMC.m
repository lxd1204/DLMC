function [label,G] = SFMC(Xs, nCluster, nAnchor)

nView = length(Xs);
for iView = 1:nView
    a = max(max(Xs{iView}));
    Xs{iView} = double(Xs{iView}./a);
end
projev = 1.5;
i = 1;
opt = [];
opt1. style = 1;
opt1. IterMax =50;
opt1. toy = 0;

[P1, alpha, label,G] = FastmultiCLR_v2(Xs, nCluster, nAnchor, opt1,10);
end