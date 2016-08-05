function index = manipulateClusters(index, manipulationSets, colorData, opts_fkmeans)

if ~isfield(manipulationSets,'removalSet') || isempty(manipulationSets.removalSet); removalSet = [];      else; removalSet = manipulationSets.removalSet; end;
if ~isfield(manipulationSets,'mergeSets')  || isempty(manipulationSets.mergeSets) ; mergeSets  = cell(0); else; mergeSets  = manipulationSets.mergeSets;  end;
if ~isfield(manipulationSets,'splitSet')   || isempty(manipulationSets.splitSet)  ; splitSet   = [];      else; splitSet   = manipulationSets.splitSet;   end;
if ~isfield(manipulationSets,'splitCount') || isempty(manipulationSets.splitCount); splitCount = 2;       else; splitCount = manipulationSets.splitCount; end;

newClusterCount                                     = numel(setdiff(unique(index), [-1 0], 'stable'));
% REMOVE
SVsToRecluster                                      = ismember(index, removalSet);
index(SVsToRecluster)                               = -1;
% MERGE
activeClusters                                      = setdiff(unique(index), [-1 0]);
for kk = 1:numel(mergeSets)
  mergeSets{kk}                                     = intersect(mergeSets{kk}, activeClusters, 'stable');
  index(ismember(index, mergeSets{kk}))             = min(mergeSets{kk});
end
% SPLIT
newClusterCount                                     = newClusterCount + 1;
for kk = 1:numel(splitSet)
  SVsToSplit                                        = find(index == splitSet(kk));
%  T                                                 = cluster(linkage(colorData(SVsToSplit, :), 'ward'), 'maxclust', 2);
  thisOpts                                          = opts_fkmeans;
  thisOpts.weight                                   = thisOpts.weight(SVsToSplit);
  T                                                 = snrAwareKmeans(splitCount, colorData(SVsToSplit, :), thisOpts);
  for mm = splitCount:-1:2
    index(SVsToSplit(T==mm))                        = newClusterCount;
    newClusterCount                                 = newClusterCount + 1;
  end
end
% CONSOLIDATE
newClusterCount                                     = numel(setdiff(unique(index), [-1 0], 'stable'));
for kk = 1:newClusterCount
  while ~any(index==kk)
    index(index>kk)                                 = index(index>kk) - 1;
  end
end


