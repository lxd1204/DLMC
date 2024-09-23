function res = ComputeECVIs(Ytrue, Ypred)
% Input:
%         Ytrue, nSmp * 1, 1 to nClass
%         Ypred, nSmp * 1, 1 to nCluster
% Output
%         Acc
%         NMI_sqrt
%         NMI_max
%         AMI
%         ARI
%         Purity
% 
% Note
% The number of clusters in Ypred could be less or larger than the number of classes in Ytrue
% 
% 

%**********************************************
% relabel Ypred into 1 to nCluster
%**********************************************
[~, ~, Ypred] = unique(Ypred); 
Cont = ComputeContingency(Ytrue, Ypred);

nSmp = length(Ytrue);


%**********************************************
% Compute Acc without nClass == nCluster
%**********************************************
M = matchpairs(Cont, 0, 'max');
acc = sum(Cont(sub2ind(size(Cont), M(:,1), M(:,2))))/nSmp;

%**********************************************
% Compute NMI without nClass == nCluster
%**********************************************
[nClass, nCluster] = size(Cont);
p_row = sum(Cont, 2)'; %  1 * nClass
p_column = sum(Cont, 1); % 1 * nCluster

% Calculating the Entropies
H_row = -(p_row / nSmp) * log2(p_row / nSmp)';
H_column = -(p_column / nSmp) * log2(p_column / nSmp)';

% Calculate the MI (unadjusted)
MI = 0;
for iClass = 1:nClass
    for iCluster = 1:nCluster
        if Cont(iClass, iCluster) > 0
            MI = MI + Cont(iClass, iCluster) * log2(Cont(iClass, iCluster) * nSmp / (p_row(iClass) * p_column(iCluster)));
        end
    end
end
MI = MI / nSmp;
nmi_sqrt = MI / max(sqrt(H_row * H_column), eps); % to avoid nCluster == 1
nmi_max = MI / max(H_row, H_column);
nmi_half = MI / ((H_row + H_column) * 0.5);
%**********************************************
% Compute ARI without nClass == nCluster
%**********************************************
nis = sum(sum(Cont, 2).^2); % Sum of squares of sums of rows
njs = sum(sum(Cont, 1).^2); % Sum of squares of sums of columns

t1 = nchoosek(nSmp, 2); % Total number of pairs of entities
t2 = sum(sum(Cont.^2)); % Sum over rows and columns of nij^2
t3 = 0.5 * (nis + njs);

% Expected index (for adjustment)
nc = (nSmp * (nSmp^2 + 1) - (nSmp + 1) * nis - (nSmp + 1) * njs + 2 * (nis * njs) / nSmp) / (2 * (nSmp - 1));

A = t1 + t2 - t3; % No. agreements

if t1 == nc
    ari = 0; % Avoid division by zero; if k=1, define Rand = 0
else
    ari = (A - nc) / (t1 - nc); % Adjusted Rand - Hubert & Arabie 1985
end

%**********************************************
% Compute Purity without nClass == nCluster
%**********************************************
n_k = max(Cont, [], 1);
purity = sum(n_k) / nSmp;

%**********************************************
% Compute AMI without nClass == nCluster
%**********************************************
% Correcting for agreement by chance
AB = p_row' * p_column;
bound = zeros(nClass, nCluster);

E3 = (AB / nSmp^2) .* log2(AB / nSmp^2);

EPLNP = zeros(nClass, nCluster);
LogNij = log2((1:min(max(p_row), max(p_column))) / nSmp);
for iClass= 1:nClass
    for iCluster = 1:nCluster
        % sumPnij = 0;
        nij = max(1, p_row(iClass) + p_column(iCluster) - nSmp);
        X = sort([nij, nSmp - p_row(iClass) - p_column(iCluster) + nij]);
        if nSmp - p_column(iCluster) > X(2)
            nom = [(p_row(iClass) - nij + 1:p_row(iClass)), (p_column(iCluster) - nij + 1:p_column(iCluster)), (X(2) + 1:nSmp - p_column(iCluster))];
            dem = [(nSmp - p_row(iClass) + 1:nSmp), (1:X(1))];
        else
            nom = [(p_row(iClass) - nij + 1:p_row(iClass)), (p_column(iCluster) - nij + 1:p_column(iCluster))];
            dem = [(nSmp - p_row(iClass) + 1:nSmp), (nSmp - p_column(iCluster) + 1:X(2)), (1:X(1))];
        end
        p0 = prod(nom ./ dem) / nSmp;
        
        sumPnij = p0;
        
        EPLNP(iClass, iCluster) = nij * LogNij(nij) * p0;
        p1 = p0 * (p_row(iClass) - nij) * (p_column(iCluster) - nij) / (nij + 1) / (nSmp - p_row(iClass) - p_column(iCluster) + nij + 1);
        
        for nij = max(1, p_row(iClass) + p_column(iCluster) - nSmp) + 1:1:min(p_row(iClass), p_column(iCluster))
            sumPnij = sumPnij + p1;
            EPLNP(iClass, iCluster) = EPLNP(iClass, iCluster) + nij * LogNij(nij) * p1;
            p1 = p1 * (p_row(iClass) - nij) * (p_column(iCluster) - nij) / (nij + 1) / (nSmp - p_row(iClass) - p_column(iCluster) + nij + 1);
        end
        CC = nSmp * (p_row(iClass) - 1) * (p_column(iCluster) - 1) / p_row(iClass) / p_column(iCluster) / (nSmp - 1) + nSmp / p_row(iClass) / p_column(iCluster);
        bound(iClass, iCluster) = p_row(iClass) * p_column(iCluster) / nSmp^2 * log2(CC);
    end
end

EMI_bound = sum(sum(bound));
% EMI_bound_2 = log2(nClass * nCluster / nSmp + (nSmp - nClass) * (nSmp - nCluster) / (nSmp * (nSmp - 1)));
EMI = sum(sum(EPLNP - E3));

ami = (MI - EMI) / (sqrt(H_row * H_column) - EMI);

% If expected mutual information is negligible, use NMI.
if abs(EMI) > EMI_bound
    fprintf('The EMI is small: EMI < %f, setting AMI = NMI\n', EMI_bound);
    ami = nmi_sqrt;
end

res = [acc; nmi_sqrt; nmi_max; nmi_half; ami; ari; purity];
end

function Cont = ComputeContingency(Ytrue, Ypred)
% 
% nSmp = length(Ytrue);
% Cont = zeros(max(Ytrue), max(Ypred));
% for iSmp = 1:length(Ytrue)
%     Cont(Ytrue(iSmp), Ypred(iSmp)) = Cont(Ytrue(iSmp), Ypred(iSmp)) + 1;
% end
Cont = accumarray([Ytrue, Ypred], 1);
end