load(mergedSvFileName);
%opts.spatialProximityRadius        = graphData.spatialNhoodRadius;
%opts.minimalComponentSize          = graphData.minSizeForPure;
%[square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, ~] = removeSmallSupervoxels(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, 1, opts);
voxelCount                          = prod(stackSize);
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0; for kk = 1:size(bbVol, 4); rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack; end; clear rawStack;


indexOrig                       = index;
svSeg = zeros(stackSize); for kk = 1:numel(svCells); svSeg(svCells{kk}) = kk; end;
groundTruth = zeros(numel(svCells), 1);
claimed = false(numel(svCells), 1);
cd /home/uygar/bb/growFlesh;
counter                               = 1;
for kk = 1:188
  fileName                            = ['/home/uygar/bb/data/dawen/nTracer/swc/voxel/PVneuron_Hippocampus Soma-1_Neuron-' num2str(kk) '_stdSWC.swc'];
  [nodes,edges,radii,nodeTypes]       = readArborComponentTrace(fileName, -1:6);
  nodes(:, 1:2)                       = nodes(:, 2:-1:1) + 1;
  ind                                 = sub2ind(stackSize, nodes(:,1), nodes(:,2), nodes(:,3));
  involvedSVs                         = setdiff(unique(svSeg(ind)), 0, 'stable');
  multiplyClaimed                     = find(claimed(involvedSVs));
  groundTruth(multiplyClaimed)        = 0;
  involvedSVs(multiplyClaimed)        = [];
  if ~isempty(involvedSVs)
    claimed(involvedSVs)              = true;
    groundTruth(involvedSVs)          = counter;
    counter                           = counter + 1;
  end
end
groundTruthOrig                       = groundTruth;
cd /home/uygar/bb;
svSizes                               = cellfun(@numel,svCells); svSizes = svSizes(:);
avgColors                             = zeros(188, size(graphData.colorsForDistal, 2));
for kk = 1:max(groundTruth)
  theseSVs                            = find(groundTruth==kk);
  avgColors(kk, :)                    = svSizes(theseSVs)' * graphData.colorsForDistal(theseSVs, :) / sum(svSizes(theseSVs));
end
cd /home/uygar/bb/quantification/clusterTools
myLinkage = elinkage(pdist(avgColors,'euclidean')); myLinkage(:,3) = myLinkage(:,3)/myLinkage(end,3);
L=zeros(stackSize); for kk=1:numel(svCells); L(svCells{kk})=kk; end;
for kk=1:graphData.opts_irbleigs.K
  tmp=false(stackSize); tmp(cell2mat(svCells(index==kk)'))=true; CC=bwconncomp(tmp);
  for mm=1:CC.NumObjects
    involvedSVs = unique(L(CC.PixelIdxList{mm}));
    if numel(CC.PixelIdxList{mm})<0.5*sum(voxelCounts(index==kk))/CC.NumObjects
      index(involvedSVs)=0;
    end
  end
end
groundTruth                           = groundTruthOrig;
groundTruth(index==0)                 = 0;
cd ~/bb/quantification/clusterTools/
hardClusterStats = cutDendrogram(myLinkage, groundTruth, index, 0.1);


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
%writeFileName = ['~/bb/results/pvhippoPerimeter_' num2str(graphData.c) '_' num2str(graphData.s) '_' num2str(graphData.colorRadiusForPure) '_' num2str(graphData.minSizeForPure) '_' num2str(graphData.maxPerim) '_'];
writeFileName = ['~/bb/results/pvhippoPerimeter_' num2str(graphData.c) '___' num2str(graphData.colorRadiusForPure) '_' num2str(graphData.minSizeForPure) '_' num2str(graphData.maxPerim) '_'];
writeFileName = [writeFileName num2str(graphData.spatialNhoodRadius) '___' num2str(graphData.maxColorRadiusForProximal) '___'];
%writeFileName = [writeFileName num2str(graphData.minEdgeCountForProximal) '_' num2str(graphData.maxAssumedSdist) '_' num2str(graphData.opts_irbleigs.K) '_ari' num2str(round(10000*max(hardClusterStats.ari))) '.jpg'];
writeFileName = [writeFileName num2str(graphData.minEdgeCountForProximal) '___' num2str(graphData.opts_irbleigs.K) '_ari' num2str(round(10000*max(hardClusterStats.ari))) '.jpg'];
imwrite(bigIm(:,:,1:3), writeFileName);
