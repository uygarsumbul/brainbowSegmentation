function [square_sAff, svMeans, svCells, voxelCounts, tmpC] = mergeConnComp3(index, square_sAff, svMeans, svCells, voxelCounts, opts, connTh, careful)

svMeansNorm                             = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
if size(svMeans,2)==3
  svMeansNormLUV                        = rgb2luv(svMeansNorm')';
else
  svMeansNormLUV(:,1:3)                 = rgb2luv(svMeansNorm(:,1:3)')';
  svMeansNormLUV(:,4)                   = svMeansNorm(:,4)*50;
end

newsvCells                              = cell(0);
binsaff                                 = (square_sAff > connTh);
tmpC                                    = zeros(1, numel(svCells));
idx                                     = 1;
for kk = 1:max(index)
  thisCluster                           = find(index==kk);
  [S, C]                                = graphconncomp(binsaff(thisCluster, thisCluster));
  for mm = 1:S
    thisConnComp                        = thisCluster(C==mm);
    ntcc                                = numel(thisConnComp);
    if careful && ntcc>1
      tmp                               = reshape(harmmean([kron(ones(ntcc,1),voxelCounts(thisConnComp)) kron(voxelCounts(thisConnComp),ones(ntcc,1))],2),ntcc,ntcc);
      [Stmp, Ctmp]                      = graphconncomp((squareform(pdist(svMeansNorm(thisConnComp,:))) < opts.maxColorDistSizeProduct ./ tmp) & binsaff(thisConnComp, thisConnComp));
      for nn = 1:Stmp
        thisConnCompTmp                 = thisConnComp(Ctmp==nn);
        newsvCells{end+1}               = cell2mat(svCells(thisConnCompTmp)');
        tmpC(thisConnCompTmp)           = idx;
        idx                             = idx + 1;
      end
    else
      newsvCells{end+1}                 = cell2mat(svCells(thisConnComp)');
      tmpC(thisConnComp)                = idx;
      idx                               = idx + 1;
    end
  end
end

S                                       = idx - 1;
[row, col, val]                         = find(square_sAff);
newRows                                 = tmpC(row);
newCols                                 = tmpC(col);
tmp                                     = find(newRows<=newCols);
newRows(tmp)                            = [];
newCols(tmp)                            = [];
val(tmp)                                = [];
newsvMeans                              = zeros(S, size(svMeans, 2));
newvoxelCounts                          = zeros(S, 1);
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
for kk = 1:S
  newvoxelCounts(kk)                    = numel(newsvCells{kk});
  thisGroup                             = find(tmpC==kk);
  newsvMeans(kk, :)                     = voxelCounts(thisGroup)' * svMeans(thisGroup, :) / sum(voxelCounts(thisGroup));
end
svMeans                                 = newsvMeans;
svCells                                 = newsvCells;
voxelCounts                             = newvoxelCounts;

binsaff                                 = (square_sAff>1/sqrt(3)-1e-5);
binsaff                                 = (sum(binsaff,2)>0);
tmp                                     = find(~binsaff & (voxelCounts<opts.disconnectedSVsizeTh));
square_sAff(tmp, :)                     = [];
square_sAff(:, tmp)                     = [];
svMeans(tmp, :)                         = [];
svCells(tmp)                            = [];
voxelCounts(tmp)                        = [];
