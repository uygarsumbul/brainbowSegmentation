function superVoxelCells = newSplitHardSuperVoxels(splitHardSVopts, superVoxelCells, bbVol)

stackSize                                                 = size(bbVol);
stackSize                                                 = stackSize(1:3);
voxelCount                                                = prod(stackSize);
detcov                                                    = zeros(1, numel(superVoxelCells));
counts                                                    = zeros(1, numel(superVoxelCells));
normalizer                                                = sqrt(sum(bbVol.^2, 4));
for dd = 1:size(bbVol, 4)
  bbVol(:,:,:,dd) = bbVol(:,:,:,dd) ./ normalizer;
end; clear normalizer;
for kk = 1:numel(superVoxelCells)
  tmp                                                     = zeros(numel(superVoxelCells{kk}), size(bbVol, 4));
  for dd = 1:size(bbVol, 4)
    tmp(:, dd)                                            = bbVol(superVoxelCells{kk} + (dd-1)*voxelCount);
  end
  detcov(kk)                                              = det(cov(tmp));
end
counts                                                    = cellfun(@numel, superVoxelCells);
detcov(counts<2)                                          = 0;
hardSuperVoxels                                           = find(detcov>splitHardSVopts.detThreshold & counts>splitHardSVopts.subdivisionSizeThreshold);
flag                                                      = true;
while ~isempty(hardSuperVoxels) & flag
  prevSVcount                                             = numel(superVoxelCells);
  newCells                                                = cell(1, numel(hardSuperVoxels));
  newSVcount                                              = 0;
  for kk = 1:numel(hardSuperVoxels)
    allVox                                                = superVoxelCells{hardSuperVoxels(kk)};
    % SUBDIVIDE BASED ON COLOR
    normalizedColorData                                   = [];
    for dd = 1:size(bbVol, 4)
      normalizedColorData                                 = [normalizedColorData bbVol(allVox + (dd-1)*voxelCount)];
    end
    T                                                     = cluster(linkage(normalizedColorData,'ward'),'maxclust',2);
    % SUBDIVIDE EACH PIECE AGAIN BASED ON CONNECTIVITY
    CC1                                                   = connectedcomponentsInLocalCoordinates(stackSize, allVox(T==1), splitHardSVopts.connectivity);
    CC2                                                   = connectedcomponentsInLocalCoordinates(stackSize, allVox(T==2), splitHardSVopts.connectivity);
    newCells{kk}                                          = [CC1.PixelIdxList CC2.PixelIdxList];
    newSVcount                                            = newSVcount + numel(CC1.PixelIdxList) + numel(CC2.PixelIdxList);
  end
  superVoxelCells(hardSuperVoxels)                        = [];
  detcov(hardSuperVoxels)                                 = [];
  counts(hardSuperVoxels)                                 = [];
  oldSvCount                                              = numel(detcov);
  detcov                                                  = [detcov zeros(1, newSVcount)];
  for kk = 1:numel(newCells)
    superVoxelCells                                       = [superVoxelCells newCells{kk}];
  end
  for kk = oldSvCount+1:oldSvCount+newSVcount
    tmp                                                   = zeros(numel(superVoxelCells{kk}), size(bbVol, 4));
    for dd = 1:size(bbVol, 4)
      tmp(:, dd)                                          = bbVol(superVoxelCells{kk} + (dd-1)*voxelCount);
    end
    detcov(kk)                                            = det(cov(tmp));
  end
  counts                                                  = [counts cellfun(@numel, superVoxelCells(oldSvCount+1:end))];
  detcov(find(counts(oldSvCount+1:end)<2)+oldSvCount)     = 0;
  hardSuperVoxels                                         = oldSvCount + find(detcov(oldSvCount+1:end)>splitHardSVopts.detThreshold & counts(oldSvCount+1:end)>splitHardSVopts.subdivisionSizeThreshold);
  flag                                                    = (prevSVcount~=numel(superVoxelCells));
end

function myCC = connectedcomponentsInLocalCoordinates(stackSize, theseVoxels, connectivity)
[xx,yy,zz]              = ind2sub(stackSize, theseVoxels);
xSub                    = min(xx)-1;
ySub                    = min(yy)-1;
zSub                    = min(zz)-1;
xx                      = xx-xSub;
yy                      = yy-ySub;
zz                      = zz-zSub;
maxxx                   = max(xx);
maxyy                   = max(yy);
maxzz                   = max(zz);
tmp                     = false(maxxx, maxyy, maxzz);
reducedIndices          = sub2ind([maxxx, maxyy, maxzz], xx, yy, zz);
tmp(reducedIndices)     = true;
myCC                    = bwconncomp(tmp, connectivity);
for nn = 1:numel(myCC.PixelIdxList)
  [xx,yy,zz]            = ind2sub(size(tmp), myCC.PixelIdxList{nn});
  xx                    = xx+xSub;
  yy                    = yy+ySub;
  zz                    = zz+zSub;
  myCC.PixelIdxList{nn} = sub2ind(stackSize, xx, yy, zz);
end
