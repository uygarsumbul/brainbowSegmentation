function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex] = mergeSmallSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, opts)

% opts.luvColorDistanceUpperBound = 10; opts.disconnectedSVsizeTh = 20;
svMeansNorm                                 = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
allTriplets                                 = nchoosek(1:size(svMeans, 2), 3);
allColors                                   = zeros(size(svMeans, 1), 3*size(allTriplets, 1));
for kk = 1:size(allTriplets, 1)
  allColors(:, 3*kk-2:3*kk)                 = rgb2luv(svMeansNorm(:, allTriplets(kk, :))')';
end
[coeff,score,latent]                        = pca(allColors);
svMeansNormLUV                              = score(:, 1:size(svMeans, 2));
voxelCounts                                 = voxelCounts(:);
Mdl                                         = KDTreeSearcher(svMeansNormLUV);
binsaff                                     = (square_sAff>0.5);
xx                                          = zeros(nnz(binsaff), 1);
yy                                          = xx;
idx                                         = 1;
cc                                          = size(svMeansNormLUV, 1);
smallSVs                                    = find(voxelCounts<opts.disconnectedSVsizeTh);
for smallSV = 1:numel(smallSVs) % VERY INEFFICIENT
  kk                                        = smallSVs(smallSV);
  [colorNeighborsCell, Dcell]               = rangesearch(Mdl, svMeansNormLUV(kk, :), opts.luvColorDistanceUpperBound);
  colorNeighbors                            = colorNeighborsCell{1}(2:end);
  D                                         = Dcell{1}(2:end);
  [colorNeighbors, ia, ~]                   = intersect(colorNeighbors, find(binsaff(kk, :)));
  D                                         = D(ia);
  [~, pos]                                  = min(D);
  bestN                                     = colorNeighbors(pos);
  maxVoxColorDist                           = max(max(svColorMaxs([kk bestN], :), [], 1) - min(svColorMins([kk bestN], :), [], 1));
  if ~isempty(pos) && maxVoxColorDist<opts.maxVoxColorDist
    yy(idx)                                 = bestN;
    xx(idx)                                 = kk;
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
