function [L, superVoxelCells] = splitHardSuperVoxels(splitHardSVopts, superVoxelCells, bbVol)

stackSize  = size(bbVol);
stackSize  = stackSize(1:3);
voxelCount = prod(stackSize);
detcov     = zeros(1, numel(superVoxelCells));
counts     = zeros(1, numel(superVoxelCells));
for kk=1:numel(superVoxelCells)
  if numel(superVoxelCells{kk})>1
    tmp        = [];
    for dd = 1:size(bbVol, 4)
      tmp      = [tmp bbVol(superVoxelCells{kk} + (dd-1)*voxelCount)];
    end
    tmp        = tmp ./ repmat(sqrt(sum(tmp.^2,2)),1,size(tmp,2));
    detcov(kk) = det(cov(tmp));
    counts(kk) = numel(superVoxelCells{kk});
  else
    detcov(kk) = 0;
    counts(kk) = 1;
  end
end
hardSuperVoxels = find(detcov>splitHardSVopts.detThreshold);
flag            = true;
while ~isempty(hardSuperVoxels) & flag
  prevSVcount           = numel(superVoxelCells);
  subdividedSuperVoxels = false(1, numel(hardSuperVoxels));
  newCells              = cell(1, numel(hardSuperVoxels));
  newSVcount            = 0;
  for kk = 1:numel(hardSuperVoxels)
    % DON'T SUBDIVIDE IF THE COMPONENT IS TINY
    if numel(superVoxelCells{hardSuperVoxels(kk)})>splitHardSVopts.subdivisionSizeThreshold
      subdividedSuperVoxels(kk) = true;
      allVox = superVoxelCells{hardSuperVoxels(kk)};
      % SUBDIVIDE BASED ON COLOR
      normalizedColorData   = [];
      for dd = 1:size(bbVol, 4)
        normalizedColorData = [normalizedColorData bbVol(allVox + (dd-1)*voxelCount)];
      end
      normalizedColorData   = normalizedColorData ./ repmat(sqrt(sum(normalizedColorData.^2,2)),1,size(normalizedColorData,2));
      T                     = cluster(linkage(normalizedColorData,'ward'),'maxclust',2);
      % SUBDIVIDE EACH PIECE AGAIN BASED ON CONNECTIVITY
      CC1 = connectedcomponentsInLocalCoordinates(stackSize, allVox(T==1), splitHardSVopts.connectivity);
      CC2 = connectedcomponentsInLocalCoordinates(stackSize, allVox(T==2), splitHardSVopts.connectivity);
      newCells{kk} = [CC1.PixelIdxList CC2.PixelIdxList];
      newSVcount   = newSVcount + numel(CC1.PixelIdxList) + numel(CC2.PixelIdxList);
    end
  end
  superVoxelCells(hardSuperVoxels(subdividedSuperVoxels)) = [];
  detcov(hardSuperVoxels(subdividedSuperVoxels))          = [];
  counts(hardSuperVoxels(subdividedSuperVoxels))          = [];
  detcov                                                  = [detcov zeros(1, newSVcount)];
  counts                                                  = [counts zeros(1, newSVcount)];
  svCountBeforeAddition                                   = numel(detcov) - newSVcount;
  for kk = 1:numel(newCells)
    superVoxelCells = [superVoxelCells newCells{kk}];
  end
  for kk = svCountBeforeAddition+1:svCountBeforeAddition+newSVcount
    if numel(superVoxelCells{kk})>1
      tmp = [];
      for dd = 1:size(bbVol, 4)
        tmp = [tmp bbVol(superVoxelCells{kk} + (dd-1)*voxelCount)];
      end
      tmp = tmp ./ repmat(sqrt(sum(tmp.^2,2)),1,size(tmp,2));
      detcov(kk) = det(cov(tmp));
      counts(kk) = numel(superVoxelCells{kk});
    else
      detcov(kk) = 0;
      counts(kk) = 1;
    end
  end
  hardSuperVoxels = find(detcov>splitHardSVopts.detThreshold);
  flag = (prevSVcount~=numel(superVoxelCells));
end

L = zeros(stackSize);
for kk = 1:numel(superVoxelCells)
  L(superVoxelCells{kk}) = kk;
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
