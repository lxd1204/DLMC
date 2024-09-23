function U = CutBP(B, nCluster, cut_type, eig_type, norm_type)

[nSmp, nAnchor] = size(B);
if nAnchor < nCluster
    error('Need more columns!');
end

if ~exist('cut_type', 'var')
    cut_type = 1;
end

if ~exist('eig_type', 'var')
    eig_type = 2;
end
if ~exist('norm_type', 'var')
    norm_type = 2;
end

switch cut_type
    case 1
        U = bcut_tcut(B, nCluster + 1);
    case 2
        U = bcut_lsc(B, nCluster + 1);
    case 3
        U = bcut_uspec(B, nCluster + 1);
end
U = real(U);
switch eig_type
    case 1
        U = U(:, 1:nCluster);
    case 2
        U = U(:, 2:nCluster + 1);
end
    
switch norm_type
    case 1
        % un-normalized
    case 2
        U = bsxfun(@rdivide, U, sqrt(max(sum(U.^2, 2), 1e-10)));
end