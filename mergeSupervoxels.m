function mergeSupervoxels(opts)

spatialDistanceCalculation.upperBound                             = 10;
load(opts.loadFilename);
svCells                                                           = superVoxelCells; clear superVoxelCells;
svMeans                                                           = superVoxelMeans; clear superVoxelMeans;
cc                                                                = numel(svCells);
[ii_sAff, ~, ss_sAff]                                             = find(sAff);
yy                                                                = ceil( (2*cc-1 - sqrt((2*cc-1)^2-8*ii_sAff))/2 );
xx                                                                = ii_sAff - cc*yy + cc + yy.*(yy+1)/2;
square_sAff                                                       = sparse(xx, yy, ss_sAff, cc, cc);
square_sAff                                                       = square_sAff + transpose(square_sAff);
clear cc; clear ii_sAff; clear ss_sAff; clear xx; clear yy; clear sAff;
opts.mergeWRTnAo.zAnisotropy                                      = opts.zAnisotropy;
voxelCounts                                                       = cellfun(@numel, svCells);
voxelCount                                                        = prod(stackSize);
oldCount                                                          = numel(svCells);
svCount                                                           = oldCount;
iterCount                                                         = 1;
dups                                                              = [];
% HEURISTICS TO MERGE SUPERVOXELS
while true
  opts.mergeSingleNeighborSuperVoxels.maxSizeForSingleNeighborSVs = quantile(voxelCounts, 0.5);
  [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~, ~, dups] = demixSupervoxels(                square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, opts.demix, 1, 1);                                          svCount(end+1) = numel(svCells);
  [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~]          = mergeWRTneighborsAndOrientations(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, 1, opts.mergeWRTnAo.normFlag, stackSize, opts.mergeWRTnAo); svCount(end+1) = numel(svCells);
  [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~]          = mergeSingleNeighborSuperVoxels(  square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, 1, opts.mergeSingleNeighborSuperVoxels);                    svCount(end+1) = numel(svCells);
  [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~]          = mergeCloseNeighborhoods(         square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, 1, opts.mergeCloseNeighborhoods);                           svCount(end+1) = numel(svCells);
  [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~]          = mergeSmallSuperVoxels(           square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, 1, opts.mergeSmallSuperVoxels);                             svCount(end+1) = numel(svCells);
  newCount = numel(svCells); if newCount==oldCount; break; end;
  oldCount = newCount; disp(svCount); iterCount = iterCount+1; dupCount = size(dups,1);
  thisFN = [opts.saveFileName '_iter' num2str(iterCount-1) '.mat'];
  save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'dupCount', 'svCount', 'svColorMins', 'svColorMaxs', '-v7.3');
end
% SPATIAL AFFINITY CALCULATION
stackSize                                                = stackSize; % for parfor
cc                                                       = numel(svCells);
boundaryVoxelsSub                                        = cell(1, cc);
parfor kk = 1:cc
  [xx,yy,zz]                                             = ind2sub(stackSize, svCells{kk});
  xSub                                                   = min(xx)-2;
  ySub                                                   = min(yy)-2;
  zSub                                                   = min(zz)-2;
  xx                                                     = xx-xSub;
  yy                                                     = yy-ySub;
  zz                                                     = zz-zSub;
  maxxx                                                  = max(xx);
  maxyy                                                  = max(yy);
  maxzz                                                  = max(zz);
  tmp                                                    = false(maxxx+1, maxyy+1, maxzz+1);
  reducedIndices                                         = sub2ind([maxxx+1, maxyy+1, maxzz+1], xx, yy, zz);
  tmp(reducedIndices)                                    = true;
  localBoundaryVoxels                                    = find(tmp & ~imerode(tmp, ones(3,3,3)));
  if ~isempty(localBoundaryVoxels)
    [xx,yy,zz]                                           = ind2sub(size(tmp), localBoundaryVoxels);
  end
  xx                                                     = xx+xSub;
  yy                                                     = yy+ySub;
  zz                                                     = zz+zSub;
  boundaryVoxelsSub{kk}                                  = [xx,yy,zz*opts.zAnisotropy];
end
sAff                                                     = calculate_sAff(cc, boundaryVoxelsSub, spatialDistanceCalculation);
[ii_sAff, ~, ss_sAff]                                    = find(sAff);
yy                                                       = ceil( (2*cc-1 - sqrt((2*cc-1)^2-8*ii_sAff))/2 );
xx                                                       = ii_sAff - cc*yy + cc + yy.*(yy+1)/2;
square_sAff                                              = sparse(xx, yy, ss_sAff, cc, cc);
square_sAff                                              = square_sAff + transpose(square_sAff);
clear ii_sAff; clear ss_sAff; clear xx; clear yy; clear sAff;
thisFN = [opts.saveFileName '_sAff' num2str(spatialDistanceCalculation.upperBound) '.mat'];
save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'dupCount', 'svCount', 'svColorMins', 'svColorMaxs', '-v7.3');
