function [square_sAff, svMeans, svCells, voxelCounts] = mergeClusters(C, noMerging, square_sAff, svMeans, svCells, voxelCounts)

% C = zeros(size(index)); C(smallSVs) = 1:numel(smallSVs); for kk=1:max(index); C(index==kk)=numel(smallSVs)+kk; end;
voxelCounts         = voxelCounts(:);
S                   = max(C);
newsvCells          = cell(1, S);
newsvMeans          = zeros(S, size(svMeans, 2));
for kk = 1:S
  thisMergedSVs     = find(C==kk);
  newsvCells{kk}    = cell2mat(svCells(thisMergedSVs)');
  newsvMeans(kk, :) = voxelCounts(thisMergedSVs)' * svMeans(thisMergedSVs, :) / sum(voxelCounts(thisMergedSVs));
end
svCells             = newsvCells;
svMeans             = newsvMeans;
voxelCounts         = cellfun(@numel, newsvCells);

[row, col, val]     = find(square_sAff);
row                 = C(row);
col                 = C(col);
noMerging           = C(noMerging);   % added line on 14.03.2016
upper               = find(row<=col); % added =    on 14.03.2016
row(upper)          = [];
col(upper)          = [];
val(upper)          = [];
ff                  = max(row)+1;
noMergingEntries    = find(ismember(row, noMerging) & ismember(col, noMerging));
mergingEntries      = setdiff(1:numel(row), noMergingEntries);
noMergingRows       = row(noMergingEntries);
noMergingCols       = col(noMergingEntries);
noMergingVals       = val(noMergingEntries);
mergingRows         = row(mergingEntries);
mergingCols         = col(mergingEntries);
mergingVals         = val(mergingEntries);
id                  = mergingRows + mergingCols*ff;
[uniqueid, ia, ~]   = unique(id);
newRows             = mergingRows(ia);
newCols             = mergingCols(ia);
newVals             = zeros(size(newRows));
parfor kk = 1:numel(uniqueid)
  newVals(kk)       = max(mergingVals(id==uniqueid(kk)));
end



square_sAff         = sparse([noMergingRows(:); newRows(:)], [noMergingCols(:); newCols(:)], [noMergingVals(:); newVals(:)], S, S);
square_sAff         = square_sAff + transpose(square_sAff);
