function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex] = mergeCloseNeighborhoods(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, opts)

svMeansNorm                                 = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
allTriplets                                 = nchoosek(1:size(svMeans, 2), 3);
allColors                                   = zeros(size(svMeans, 1), 3*size(allTriplets, 1));
for kk = 1:size(allTriplets, 1)
  allColors(:, 3*kk-2:3*kk)                 = rgb2luv(svMeansNorm(:, allTriplets(kk, :))')';
end
[coeff,score,latent]                        = pca(allColors);
svMeansNormLUV                              = score(:, 1:size(svMeans, 2));
binsaff                                     = (square_sAff>1/(sqrt(3)+1e-5));
xx                                          = zeros(nnz(binsaff), 1);
yy                                          = xx;
idx                                         = 1;
cc                                          = size(svMeans, 1);
voxelCounts                                 = voxelCounts(:);

for kk = 1:cc
  neighbors = find(binsaff(kk, :));
  if numel(neighbors)>1 & max(pdist(svMeansNormLUV(neighbors, :)))<opts.maxDistNormLUV & max(svColorMaxs(neighbors, :), [], 1)-min(svColorMins(neighbors, :), [] , 1)<opts.maxVoxColorDist
    xx(idx:idx+numel(neighbors)-1)          = kk;
    yy(idx:idx+numel(neighbors)-1)          = neighbors;
    idx                                     = idx + numel(neighbors);
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

