function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex] = removeSmallSupervoxels(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, opts)

binsaff                                   = (square_sAff>=1/opts.spatialProximityRadius);
[S, C]                                    = graphconncomp(binsaff);
toRemove                                  = [];
for kk = 1:numel(svCells)
  thisComponent                           = find(C==kk);
  if sum(voxelCounts(thisComponent))<opts.minimalComponentSize
    toRemove                              = [toRemove thisComponent];
  end
end
square_sAff(toRemove, :)                  = [];
square_sAff(:, toRemove)                  = [];
svMeans(toRemove, :)                      = [];
svCells(toRemove)                         = [];
svColorMins(toRemove, :)                  = [];
svColorMaxs(toRemove, :)                  = [];
voxelCounts(toRemove)                     = [];
toRemove2                                 = ismember(origIndex, toRemove);
origIndex(toRemove2)                      = [];
mapDown                                   = zeros(1, max(origIndex));
uniqueIndices                             = unique(origIndex); % sorts by default
mapDown(uniqueIndices)                    = 1:numel(uniqueIndices);
origIndex                                 = mapDown(origIndex);
