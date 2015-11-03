function index = snrAwareRetainLargeClumps(square_sAff, index, minClumpCount, minSVcountFraction, clustersToOperate)

oversegmentClusterCount                             = max(index);
if nargin<5
  clustersToOperate                                 = 1:oversegmentClusterCount;
end
toRemove                                            = [];
for mm = 1:numel(clustersToOperate)
  kk                                                = clustersToOperate(mm);
  thisSubset                                        = find(index==kk);
  localGraph                                        = square_sAff(thisSubset, thisSubset);
  localGraph                                        = localGraph .* (localGraph>1/(sqrt(3)+1e-6));
  [S, C]                                            = graphconncomp(localGraph);
  tt                                                = hist(C,1:S);
  smallClumpIDs                                     = find(tt<minClumpCount);
  toRemove                                          = [toRemove; thisSubset(ismember(C, smallClumpIDs))];
end
index(toRemove)                                     = 0;
for kk = 1:oversegmentClusterCount
  if nnz(index==kk)<minSVcountFraction*nnz(index>0)
    index(index==kk)                                = 0;
  end
end
% CONSOLIDATE
newClusterCount                                     = numel(setdiff(unique(index), [-1 0], 'stable'));
for kk = 1:newClusterCount
  while ~any(index==kk)
    index(index>kk)                                 = index(index>kk) - 1;
  end
end
