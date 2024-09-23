function paramCell = EOMSC_build_param(anchor_sizes, d_rates)

if ~exist('anchor_sizes', 'var')
    anchor_sizes = [1:7];
end

if ~exist('d_rates', 'var')
    d_rates = [1.6];
end

nParam = length(anchor_sizes) * length(d_rates);
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(anchor_sizes)
    for i2 = 1:length(d_rates)
        param = [];
        param.anchors = anchor_sizes(i1);
        param.d_rate = d_rates(i2);
        idx = idx + 1;
        paramCell{idx,1} = param;
    end
end
end