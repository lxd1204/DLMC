function [F,label] = EMKMC_Bs(Bs, nCluster, gamma)
if ~exist('gamma', 'var')
    % 1<gamma<2
    gamma = 1.6;
end
maxIter = 10;
nView = length(Bs);
nSmp = size(Bs{1}, 1);
alpha = ones(nView,1)/nView;
F = initialize(nSmp,nCluster);
Gs = cell(1, nView);
for v = 1:nView
    Gs{v} = rand(size(Bs{v}, 2),nCluster);
end
for Iter = 1:maxIter
    %*******************************
    % Update F
    %*******************************
    Gtemp = 0;
    for v = 1:nView
        Gtemp=Gtemp + (alpha(v)^gamma)*Bs{v}*Gs{v};
    end
    [AA,~,CC] = svd(Gtemp','econ');
    F = (AA*CC')';
    %*******************************
    % Update G{v}
    %*******************************
    %Ftemp = F*pinv(F'*F);
    for v = 1:nView
        Gs{v} = Bs{v}' * F;
    end
    %*******************************
    % Update \alpha
    %*******************************
    Wtemp = zeros(nView, 1);
    for v = 1:nView
        Ev = Bs{v}-F*Gs{v}';
        ev = sum(sum(Ev.^2));
        r = 1/(1-gamma);
        Wtemp(v) = (gamma*ev)^r;
    end
    alpha = Wtemp./(sum(Wtemp,2));
end
[~,label] = max(F,[],2);