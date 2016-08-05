function [square_sAff, svMeans, svCells, voxelCounts, origIndex] = mergeSingleNeighborSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, origIndex, opts)

binsaff                                     = (square_sAff>1/(sqrt(3)+1e-5));
xx                                          = zeros(nnz(binsaff), 1);
yy                                          = xx;
idx                                         = 1;
cc                                          = size(svMeans, 1);
singleNeighborSVs                           = find(sum(binsaff, 2)==1 & voxelCounts<=opts.maxSizeForSingleNeighborSVs);
for singleNeighborSV = 1:numel(singleNeighborSVs)
  kk                                        = singleNeighborSVs(singleNeighborSV);
  xx(idx)                                   = kk;
  yy(idx)                                   = find(binsaff(kk, :));
  idx                                       = idx + 1;
end
xx(idx:end)                                 = [];
yy(idx:end)                                 = [];
binsaff                                     = sparse(xx, yy, 1, cc, cc);
[S,C]                                       = graphconncomp(binsaff, 'Weak', true);

newsvCells                                  = cell(1, S);
newsvMeans                                  = zeros(S, size(svMeans, 2));
for kk = 1:S
  thisConnComp                              = find(C==kk);
  newsvCells{kk}                            = cell2mat(svCells(thisConnComp)');
  newsvMeans(kk, :)                         = voxelCounts(thisConnComp)' * svMeans(thisConnComp, :) / sum(voxelCounts(thisConnComp));
end
newvoxelCounts                              = cellfun(@numel, newsvCells);

[row, col, val]                         = find(square_sAff);
newRows                                 = C(row);
newCols                                 = C(col);
tmp                                     = find(newRows<=newCols);
newRows(tmp)                            = [];
newCols(tmp)                            = [];
val(tmp)                                = [];
rrr                                     = zeros(size(row));
ccc                                     = rrr;
vvv                                     = rrr;
hhh                                     = rrr;
thisPos                                 = 1;
clear row; clear col;
for kk = 1:numel(newRows)
  thisHash                              = newRows(kk)*(S+1) + newCols(kk);
  pos                                   = find(thisHash==hhh);
  if isempty(pos)
    rrr(thisPos)                        = newRows(kk);
    ccc(thisPos)                        = newCols(kk);
    vvv(thisPos)                        = val(kk);
    hhh(thisPos)                        = thisHash;
    thisPos                             = thisPos + 1;
  else
    vvv(pos)                            = max(vvv(pos), val(kk));
  end
end
clear val;
unused                                  = find(rrr==0, 1);
rrr(unused:end)                         = [];
ccc(unused:end)                         = [];
vvv(unused:end)                         = [];
clear hhh;
square_sAff                             = sparse(rrr, ccc, vvv, S, S);
square_sAff                             = square_sAff + transpose(square_sAff);
svMeans                                 = newsvMeans;
svCells                                 = newsvCells;
voxelCounts                             = newvoxelCounts';

origIndex                               = C(origIndex);

%binsaff                                 = (square_sAff>1/sqrt(3)-1e-5);
%binsaff                                 = (sum(binsaff,2)>0);
%tmp                                     = find(~binsaff & (voxelCounts<opts.disconnectedSVsizeTh));
%square_sAff(tmp, :)                     = [];
%square_sAff(:, tmp)                     = [];
%svMeans(tmp, :)                         = [];
%svCells(tmp)                            = [];
%voxelCounts(tmp)                        = [];
