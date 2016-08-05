function [square_sAff, svMeans, svCells, voxelCounts, tmpC] = mergeSubgraphs2(square_sAff, svMeans, voxelCounts, disconnectedSVsizeTh, subgraph_svCells, subgraph_indices, nodePartition)

% nodePartition: cell variable holding the supervoxels given to individual subsets of the node partition
% subgraph_indices: cell variable holding the index variables for individual subsets of the node parition (returned by mergeConnectedComponentsOfSubgraph)
% subgraph_svCells: cell variable holding the supervoxel cells for individual subsets of the node parition (returned by mergeConnectedComponentsOfSubgraph)

oldsvCount                                  = size(square_sAff, 1);
newsvCount                                  = sum(cellfun(@numel, subgraph_svCells));
newsvCells                                  = cell(1, newsvCount);
tmpC                                        = zeros(1, oldsvCount);
offset                                      = 0;
for kk = 1:numel(subgraph_indices)
  tmpC(nodePartition{kk})                   = offset + subgraph_indices{kk};
  newsvCells(offset+1:offset+numel(subgraph_svCells{kk})) = subgraph_svCells{kk};
  offset                                    = offset + numel(subgraph_svCells{kk});
end

[row, col, val]                             = find(square_sAff);
newRows                                     = tmpC(row);
newCols                                     = tmpC(col);
tmp                                         = find(newRows<=newCols);
newRows(tmp)                                = [];
newCols(tmp)                                = [];
val(tmp)                                    = [];
clear row; clear col;
newsvMeans                                  = zeros(newsvCount, size(svMeans, 2));
newvoxelCounts                              = zeros(newsvCount, 1);
rrr                                         = zeros(size(newRows));
ccc                                         = rrr;
vvv                                         = rrr;
hhh                                         = rrr;
thisPos                                     = 1;
for kk = 1:numel(newRows)
  thisHash                                  = newRows(kk)*(newsvCount+1) + newCols(kk);
  pos                                       = find(thisHash==hhh);
  if isempty(pos)
    rrr(thisPos)                            = newRows(kk);
    ccc(thisPos)                            = newCols(kk);
    vvv(thisPos)                            = val(kk);
    hhh(thisPos)                            = thisHash;
    thisPos                                 = thisPos + 1;
  else
    vvv(pos)                                = max(vvv(pos), val(kk));
  end
end
clear val;
unused                                      = find(rrr==0, 1);
rrr(unused:end)                             = [];
ccc(unused:end)                             = [];
vvv(unused:end)                             = [];
clear hhh;
square_sAff                                 = sparse(rrr, ccc, vvv, newsvCount, newsvCount);
square_sAff                                 = square_sAff + transpose(square_sAff);
for kk = 1:newsvCount
  newvoxelCounts(kk)                        = numel(newsvCells{kk});
  thisGroup                                 = find(tmpC==kk);
  newsvMeans(kk, :)                         = voxelCounts(thisGroup)' * svMeans(thisGroup, :) / sum(voxelCounts(thisGroup));
end
svMeans                                     = newsvMeans;
svCells                                     = newsvCells;
voxelCounts                                 = newvoxelCounts;

binsaff                                     = (square_sAff>1/sqrt(3)-1e-5);
binsaff                                     = (sum(binsaff,2)>0);
tmp                                         = find(~binsaff & (voxelCounts<disconnectedSVsizeTh));
square_sAff(tmp, :)                         = [];
square_sAff(:, tmp)                         = [];
svMeans(tmp, :)                             = [];
svCells(tmp)                                = [];
voxelCounts(tmp)                            = [];
