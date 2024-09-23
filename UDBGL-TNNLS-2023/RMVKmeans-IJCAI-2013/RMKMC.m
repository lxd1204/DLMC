function [label, centerCell, objEnd, objHistory] = RMKMC(Xcell, gamma, nCluster, maxIter1, maxIter2, nRepeat)
% Input:
%     Xlow : nSmp * nFea
%     Xup : nSmp * nFea
%     nCluster: the number of cluster
%     maxIter: the number of iterations, 100 default
%     nRepeat: the number of random initialization
%
% Output:
%     label: nSmp * 1
%     center: nCluster * nFea
%     objEnd: the final objValue
%
%
% written by Liang Du (duliang@sxu.edu.cn)
% 2021/03/8


if ~exist('maxIter1', 'var')
    maxIter1 = 100;
end

if ~exist('maxIter2', 'var')
    maxIter2 = 10;
end

if ~exist('nRepeat', 'var')
    nRepeat = 1;
end

nSmp = size(Xcell{1}, 1);
nView = length(Xcell);
f1 = cell(1, nView);
for i1 = 1:nView
    f1{i1} = sum(Xcell{i1}.^2,2);
end

bestlabel = [];
for iRepeat = 1:nRepeat
    
    objHistory = [];
    
    %****************************************
    % init Kmeans center and partition
    %****************************************
    label = randi(nCluster, nSmp, 1);
    E = sparse((1:nSmp)', label, 1, nSmp, nCluster, nSmp);
    tmp = bsxfun(@times, E, 1./(sum(E,1)+eps));
    tmp = full(tmp);
    centerCell = cell(1, nView);
    for i1 = 1:nView
        centerCell{i1} = tmp' * Xcell{i1};
    end
    
    a = ones(nView,1)/nView;
    %****************************************
    % init scale parameter
    %****************************************
    e = mvKmeansL2Errors(Xcell, centerCell, label);
    w = cell(1,nView);
    for i1 = 1:nView
        w{i1} = 1 ./ (2 * sqrt(e{i1}) + eps);
    end
    
    lastLabel = zeros(nSmp,1);
    for iter1 = 1:maxIter1
        
        %****************************************
        % Update kmeans center and partition
        %****************************************
        iter2 = 0;
        while any(label ~= lastLabel) && iter2 < maxIter2
            lastLabel = label;
            
            %****************************************
            % Update partition
            %****************************************
            D = zeros(nSmp, nCluster);
            for i1 = 1:nView
                t1 = Xcell{i1} * centerCell{i1}';
                t2 = sum(centerCell{i1}.^2, 2);
                t3 = bsxfun(@plus, f1{i1}, t2');
                Dtmp = t3 - 2 * t1;
                Dtmp = bsxfun(@times, Dtmp, w{i1} * a(i1).^gamma);
                Dtmp (Dtmp < 0) = 0;
                D = D + Dtmp;
            end
            [val, label] = min(D,[],2); % assign samples to the nearest centers
            
            ll = unique(label);
            if length(ll) < nCluster
                missCluster = 1 : nCluster;
                missCluster(ll) = [];
                missNum = length(missCluster);
                [~, idx] = sort(val,1,'descend');
                label(idx(1:missNum)) = missCluster;
            end
            
            %****************************************
            % Update center for weighted kmeans
            %****************************************
            for i1 = 1:nView
                E = sparse(1:nSmp, label, w{i1} * a(i1).^gamma, nSmp, nCluster, nSmp);
                E = bsxfun(@times, E, 1./(sum(E,1) + eps));
                centerCell{i1} = full(E)' * Xcell{i1};
            end
            
            iter2 = iter2 + 1;
            
            if nargout > 3
                e = mvKmeansL2Errors(Xcell, centerCell, label);
                obj = 0;
                for i1 = 1:nView
                    obj = obj + a(i1) * sum(sqrt(e{i1}));
                end
                objHistory = [objHistory; obj];
            end
        end
        
        
        %****************************************
        % update scale parameter
        %****************************************
        e = mvKmeansL2Errors(Xcell, centerCell, label);
        
        h = zeros(nView, 1);
        for i1 = 1:nView
            h(i1) = sum(w{i1} .* e{i1});
        end
        
        a = (gamma * h).^ (1/(1-gamma));
        a = (a + eps)/max(sum(a), eps);
        
        w = cell(1,nView);
        for i1 = 1:nView
            w{i1} = 1 ./ (2 * sqrt(e{i1}) + eps);
        end
        
    end
    
    
    if isempty(bestlabel)
        bestlabel = label;
        bestcenterCell = centerCell;
        bestobjHistory = objHistory;
        
        e = mvKmeansL2Errors(Xcell, centerCell, label);
        obj = 0;
        for i1 = 1:nView
            obj = obj + a(i1) * sum(sqrt(e{i1}));
        end
        
        bestobjEnd = obj;
    else
        e = mvKmeansL2Errors(Xcell, centerCell, label);
        obj = 0;
        for i1 = 1:nView
            obj = obj + a(i1) * sum(sqrt(e{i1}));
        end
        objEnd = obj;
        if objEnd < bestobjEnd
            bestlabel = label;
            bestcenterCell = centerCell;
            bestobjEnd = objEnd;
            bestobjHistory = objHistory;
        end
    end
end
label = bestlabel;
centerCell = bestcenterCell;
objEnd = bestobjEnd;
objHistory = bestobjHistory;
end

function e = mvKmeansL2Errors(Xcell, centerCell, label)
nSmp = size(Xcell{1}, 1);
nView = length(Xcell);
f1 = cell(1, nView);
f2 = cell(1,nView);
f3 = cell(1, nView);
e =  cell(1, nView);
for i1 = 1:nView
    f1{i1} = sum(Xcell{i1}.^2,2);
    f2{i1} = sum(centerCell{i1}.^2,2);
    tmp = Xcell{i1} * centerCell{i1}';
    lidx_tmp = sub2ind(size(tmp), (1:nSmp)', label);
    f3{i1} = tmp(lidx_tmp);
    tmp2 = f2{i1};
    e{i1} = f1{i1} + tmp2(label) - 2 * f3{i1};
    e{i1} (e{i1} < 0) = 0;
end
end