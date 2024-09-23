function paramCell = SFMC_build_param(anchors)
if ~exist('anchors', 'var')
    anchors = 1000;
end

nParam = length(anchors) ;
paramCell = cell(nParam, 1);
idx = 0;
for i1 = 1:length(anchors)
    param = [];
    param.nAnchor = anchors(i1);
    idx = idx + 1;
    paramCell{idx,1} = param;
end
end