function [square_sAff, svMeans, svCells, voxelCounts] = mergeConnectedComponentsWithinIndividualClusters(index, square_sAff, svMeans, svCells, voxelCounts, opts, nhood)

if nhood == 6;  connectivityThreshold = 0.9; end;
if nhood == 18; connectivityThreshold = 0.7; end;
if nhood == 26; connectivityThreshold = 0.51; end;

newsvCells                                   = cell(0);
binsaff                                      = (square_sAff > connectivityThreshold);
tmpC                                         = zeros(1, numel(svCells));
idx                                          = 1;
for kk = 1:max(index)
  thisCluster                                = find(index==kk);
  [S, C]                                     = graphconncomp(binsaff(thisCluster, thisCluster));
  for mm = 1:S
    thisConnComp                             = thisCluster(find(C==mm));
    tmpC(thisConnComp)                       = idx;
    idx                                      = idx + 1;
    newsvCells{end+1}                        = [];
    for nn = 1:numel(thisConnComp)
      newsvCells{end}                        = [newsvCells{end}; svCells{thisConnComp(nn)}];
    end
  end
end

S                                            = idx - 1;
allGroups                                    = cell(S, 1);
for kk = 1:S
  allGroups{kk}                              = find(tmpC==kk);
end
newvoxelCounts                               = zeros(S, 1);
[row, col, val]                              = find(square_sAff);
newRows                                      = cell(numel(row)/2, 1);
newCols                                      = newRows;
newVals                                      = newRows;
newsvMeans                                   = zeros(S, size(svMeans, 2));
newvoxelCounts                               = zeros(S, 1);
entryCount_lowerHalf                         = 0;
parfor kk = 1:S
  neighborGroups                             = unique(tmpC(unique(col(ismember(row, allGroups{kk})))));
  neighborGroups(neighborGroups<=kk)         = [];
  newRows{kk}                                = kk * ones(numel(neighborGroups), 1);
  newCols{kk}                                = neighborGroups;
  for mm = 1:numel(neighborGroups)
    newVals{kk}(end+1)                       = max(max(square_sAff(allGroups{kk}, allGroups{neighborGroups(mm)})));
  end
  newvoxelCounts(kk)                         = numel(newsvCells{kk});
  newsvMeans(kk, :)                          = voxelCounts(allGroups{kk})' * svMeans(allGroups{kk}, :) / sum(voxelCounts(allGroups{kk}));
  entryCount_lowerHalf                       = entryCount_lowerHalf + numel(neighborGroups);
end
clear row; clear col; clear val;
thisCount = 0; for kk = 1:S; thisCount = thisCount + numel(newRows{kk}); end;
vecRows                                     = zeros(thisCount, 1);
vecCols                                     = zeros(1, thisCount);
vecVals                                     = zeros(1, thisCount);
idx                                         = 1;
for kk = 1:S
  vecRows(idx:idx+numel(newRows{kk})-1)     = newRows{kk};
  vecCols(idx:idx+numel(newRows{kk})-1)     = newCols{kk};
  vecVals(idx:idx+numel(newRows{kk})-1)     = newVals{kk};
  idx                                       = idx + numel(newRows{kk});
end
clear newRows; clear newCols; clear newVals; 
square_sAff                                 = sparse(vecRows, vecCols, vecVals, S, S);
square_sAff                                 = square_sAff + square_sAff';

svMeans                                     = newsvMeans;
svCells                                     = newsvCells;
voxelCounts                                 = newvoxelCounts;

binsaff                                     = (square_sAff>1/sqrt(3)-1e-5);
binsaff                                     = (sum(binsaff,2)>0);
tmp                                         = find(~binsaff & (voxelCounts<opts.disconnectedSVsizeTh));
square_sAff(tmp, :)                         = [];
square_sAff(:, tmp)                         = [];
svMeans(tmp, :)                             = [];
svCells(tmp)                                = [];
voxelCounts(tmp)                            = [];
