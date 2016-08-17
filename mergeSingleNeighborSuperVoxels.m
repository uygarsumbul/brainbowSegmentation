function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex] = mergeSingleNeighborSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, opts)

binsaff                                     = (square_sAff>1/(sqrt(3)+1e-5));
xx                                          = zeros(nnz(binsaff), 1);
yy                                          = xx;
idx                                         = 1;
cc                                          = size(svMeans, 1);
voxelCounts                                 = voxelCounts(:);
singleNeighborSVs                           = find(sum(binsaff, 2)==1 & voxelCounts<=opts.maxSizeForSingleNeighborSVs);
for singleNeighborSV = 1:numel(singleNeighborSVs)
  kk                                        = singleNeighborSVs(singleNeighborSV);
  neighbor                                  = find(binsaff(kk, :));
  if max(max(svColorMaxs([kk neighbor], :), [], 1)-min(svColorMins([kk neighbor], :), [], 1))<opts.maxVoxColorDist
    xx(idx)                                 = kk;
    yy(idx)                                 = neighbor;
    idx                                     = idx + 1;
  end
end
xx(idx:end)                                 = [];
yy(idx:end)                                 = [];
binsaff                                     = sparse(xx, yy, 1, cc, cc);
[S,C]                                       = graphconncomp(binsaff, 'Weak', true);

newsvCells              = cell(1, S);
newsvMeans              = zeros(S, size(svMeans, 2));
newsvColorMins          = zeros(S, size(svMeans, 2));
newsvColorMaxs          = zeros(S, size(svMeans, 2));
for kk = 1:S
  thisConnComp          = find(C==kk);
  newsvCells{kk}        = cell2mat(svCells(thisConnComp)');
  newsvMeans(kk, :)     = voxelCounts(thisConnComp)' * svMeans(thisConnComp, :) / sum(voxelCounts(thisConnComp));
  newsvColorMins(kk, :) = min(svColorMins(thisConnComp, :), [], 1);
  newsvColorMaxs(kk, :) = max(svColorMaxs(thisConnComp, :), [], 1);
end
svCells                 = newsvCells;
svMeans                 = newsvMeans;
svColorMins             = newsvColorMins;
svColorMaxs             = newsvColorMaxs;
voxelCounts             = cellfun(@numel, newsvCells);

[row, col, val]         = find(square_sAff);
row                     = C(row);
col                     = C(col);
upper                   = find(row<=col);
row(upper)              = [];
col(upper)              = [];
val(upper)              = [];

ff                      = max(row)+1;
id                      = row + col*ff;
[uniqueid, ia, ~]       = unique(id);
newRows                 = row(ia);
newCols                 = col(ia);
newVals                 = zeros(size(newRows));
parfor kk = 1:numel(uniqueid)
  newVals(kk)           = max(val(id==uniqueid(kk)));
end
square_sAff             = sparse(newRows, newCols, newVals, S, S);
square_sAff             = square_sAff + transpose(square_sAff);

origIndex               = C(origIndex);

