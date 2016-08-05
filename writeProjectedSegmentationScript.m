clusterCount                                 = max(index);
segmentation = zeros(stackSize); for kka=1:clusterCount; thisCluster = find(index==kka); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = kka; end; end;
voxelCount                                   = numel(segmentation);
xTileCount                                   = round(sqrt((clusterCount+1)/3) * 1) +2;
yTileCount                                   = ceil((clusterCount+1)/xTileCount);
compx                                        = size(segmentation,1)+size(segmentation,3);
compy                                        = size(segmentation,2)+size(segmentation,3);
bigIm                                        = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, size(bbVol,4));
xTile                                        = 1;
yTile                                        = 1;
bigIm(1:size(segmentation,1), 1:size(segmentation,2), :) = squeeze(max(bbVol, [], 3));

for kk=1:clusterCount
  [xTile, yTile]                             = ind2sub([xTileCount yTileCount], kk+1);
  tmp2                                       = find(segmentation==kk);
  for mm = 1:3
    tmp1                                     = zeros(size(segmentation));
    tmp1(tmp2)                               = bbVol(tmp2+(mm-1)*voxelCount);
    bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, mm) = max(tmp1, [], 3);
  end
end
writeFileName = ['~/bb/results/pvhippo_' num2str(graphData.c) '_' num2str(graphData.s) '_' num2str(graphData.colorRadiusForPure) '_' num2str(graphData.minSizeForPure) '_' num2str(graphData.maxDetCovForPure) '_'];
writeFileName = [writeFileName num2str(graphData.spatialNhoodRadius) '_' num2str(graphData.maxColorRadiusForPure) '_' num2str(graphData.maxColorRadiusForProximal) '_' num2str(graphData.minEdgeCountForPure) '_'];
writeFileName = [writeFileName num2str(graphData.minEdgeCountForProximal) '_' num2str(graphData.maxAssumedSdist) '_' num2str(graphData.opts_irbleigs.K) '.jpg'];
imwrite(bigIm(:,:,1:3), writeFileName);

