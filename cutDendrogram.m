function hardClusterStats = cutDendrogram(myLinkage, groundTruth, index)
% myLinkage         : linkage of ground truth traces (based on their colors)
% groundTruth       : supervoxel-based ground truth trace labels (0 if label is ambiguous)
% index             : clustering of supervoxels

% normalize linkage values
myLinkage(:,3) = myLinkage(:,3)/myLinkage(end,3);
% cutting levels for the linkage
threshold = myLinkage(:,3)-eps;
randIndices = zeros(size(threshold)); adjustedRandIndices = zeros(size(threshold)); typeConfs = zeros(size(threshold)); widths = zeros(size(threshold));
rowSplits = zeros(size(threshold)); columnSplits = zeros(size(threshold));

unambiguousPos                             = find(groundTruth);
index                                      = index(unambiguousPos);
allGT                                      = ones(numel(unambiguousPos), numel(threshold));
groundTruth                                = groundTruth(unambiguousPos);
for kk = 1:numel(threshold)
  % cut the dendrogram at threshold(kk)
  T                                        = cluster(myLinkage, 'cutoff', threshold(kk), 'criterion', 'distance');
  thisGT                                   = T(groundTruth);
  allGT(:,kk)                              = thisGT;
  [AR, RI, MI, HI, rS, cS, typeConfusions] = reportConfusionsAndRI(thisGT, index);
  randIndices(kk) = RI; adjustedRandIndices(kk) = AR; rowSplits(kk) = rS; columnSplits(kk) = cS; typeConfs(kk) = typeConfusions;
end

% calculate the width of the range of cutting levels yielding the exact same clustering, for each cutting level
for kk = 1:numel(threshold)
  width=1;while (kk-width>0);[AR,RI1,MI,HI,~,~,~]=reportConfusionsAndRI(allGT(:,kk-width),allGT(:,kk));
  if RI1==1;width=width+1;else;break;end;end;
  tmpwidth=0;while (kk+tmpwidth+1<=numel(threshold));[AR,RI2,MI,HI,~,~,~]=reportConfusionsAndRI(allGT(:,kk+tmpwidth+1),allGT(:,kk));
  if RI2==1;tmpwidth=tmpwidth+1;else;break;end;end;
  widths(kk) = width+tmpwidth;
end

[mini0, pos0] = min(typeConfs); allMin = find(typeConfs==mini0); % among minimum type-confusions ...
[~, postemp] = max(randIndices(allMin)); pos0=allMin(postemp); th0=threshold(pos0); width=widths(pos0);
myRowSplits = rowSplits(pos0); myColumnSplits = columnSplits(pos0);
T = cluster(myLinkage,'cutoff',th0,'criterion','distance');

hardClusterStats.minTypeConfs = mini0;
hardClusterStats.structuralSplits = myRowSplits;
hardClusterStats.geneticSplits = myColumnSplits;
hardClusterStats.minCutForMinTypeConfs = th0;
hardClusterStats.cutWidthForMinTypeConfs = width*cutStep;
hardClusterStats.ClusterCountAtMinCut = numel(unique(T));
hardClusterStats.hardClusters = T;
hardClusterStats.riAtMinCut = randIndices(pos0);
hardClusterStats.ariAtMinCut = adjustedRandIndices(pos0);
hardClusterStats.ri = randIndices;
hardClusterStats.ari = adjustedRandIndices;
hardClusterStats.allTypeConfs = typeConfs;
hardClusterStats.allCutWidths = widths*cutStep;
