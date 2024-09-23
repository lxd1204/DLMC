function paramCell = CutBP_build_param(cut_types, eig_types, norm_types)

if ~exist('cut_types', 'var')
    cut_types = [1, 2, 3];
end

if ~exist('eig_types', 'var')
    eig_types = [1, 2];
end

if ~exist('norm_types', 'var')
    norm_types = [1, 2];
end


nParam = length(cut_types) * length(eig_types) * length(norm_types);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(cut_types)
    for i2 = 1:length(eig_types)
        for i3 = 1:length(norm_types)
            param = [];
            param.cut_type = cut_types(i1);
            param.eig_type = eig_types(i2);
            param.norm_type = norm_types(i3);
            idx = idx + 1;
            paramCell{idx,1} = param;
        end
    end
end