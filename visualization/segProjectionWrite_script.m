clusterCount                                 = max(index);
segmentation = zeros(stackSize); for kka=1:clusterCount; thisCluster = find(index==kka); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = kka; end; end;
voxelCount                                   = numel(segmentation);
xTileCount                                   = ceil(sqrt(clusterCount/40) * 5);
yTileCount                                   = ceil(clusterCount/xTileCount);
compx                                        = size(segmentation,1)+size(segmentation,3);
compy                                        = size(segmentation,2)+size(segmentation,3);
bigIm                                        = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, 3);
for kk=1:clusterCount
  [xTile, yTile]                             = ind2sub([xTileCount yTileCount], kk);
  tmp2                                       = find(segmentation==kk);
  for mm = 1:3
    tmp1                                     = zeros(size(segmentation));
    tmp1(tmp2)                               = bbVol(tmp2+(mm-1)*voxelCount);
    bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, mm) = max(tmp1, [], 3);
  end
end
cd /home/uygar/bb/results/
imwrite(bigIm, 'thisBigIm.jpg');
cd /home/uygar/bb/
