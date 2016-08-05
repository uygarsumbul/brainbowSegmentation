function [square_sAff, svMeans, svCells, voxelCounts, origIndex] = mergeSmallSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, origIndex, opts)

% opts.luvColorDistanceUpperBound = 10; opts.disconnectedSVsizeTh = 20;
svMeansNorm                                  = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
if size(svMeans,2)==3
  svMeansNormLUV                             = rgb2luv(svMeansNorm')';
else
  svMeansNormLUV(:,1:3)                      = rgb2luv(svMeansNorm(:,1:3)')';
  svMeansNormLUV(:,4)                        = svMeansNorm(:,4)*50;
end

Mdl                                         = KDTreeSearcher(svMeansNormLUV);
binsaff                                     = (square_sAff>0.99);
xx                                          = zeros(nnz(binsaff), 1);
yy                                          = xx;
idx                                         = 1;
cc                                          = size(svMeansNormLUV, 1);
smallSVs                                    = find(voxelCounts<opts.disconnectedSVsizeTh);
for smallSV = 1:numel(smallSVs)
  kk                                        = smallSVs(smallSV);
  [colorNeighborsCell, Dcell]               = rangesearch(Mdl, svMeansNormLUV(kk, :), opts.luvColorDistanceUpperBound);
  colorNeighbors                            = colorNeighborsCell{1}(2:end);
  D                                         = Dcell{1}(2:end);
  [colorNeighbors, ia, ~]                   = intersect(colorNeighbors, find(binsaff(kk, :)));
  D                                         = D(ia);
  [~, pos]                                  = min(D);
  if ~isempty(pos)
    yy(idx)                                 = colorNeighbors(pos);
    xx(idx)                                 = kk;
    idx                                     = idx + 1;
  end
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
