function mergeSupervoxels(opts)

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
oldCount                                                          = numel(svCells);
counts_neighborsAndOrientations                                   = oldCount;
% STAGE 1 ---- HEURISTICS TO MERGE SUPERVOXELS
while true
  opts.mergeSingleNeighborSuperVoxels.maxSizeForSingleNeighborSVs = quantile(voxelCounts, 0.5);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = newMergeWRTneighborsAndOrientations(square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeWRTnAo.normFlag, stackSize, opts.mergeWRTnAo);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeSmallSuperVoxels(              square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeSmallSuperVoxels);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeSingleNeighborSuperVoxels(     square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeSingleNeighborSuperVoxels);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeCloseNeighborhoods(            square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeCloseNeighborhoods);
  newCount = numel(svCells); if newCount==oldCount; break; end;
  oldCount = newCount; counts_neighborsAndOrientations(end+1) = oldCount; disp(counts_neighborsAndOrientations);
  thisFN = [opts.saveFileName 'stage1_iter' num2str(numel(counts_neighborsAndOrientations)-1) '.mat'];
  save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'counts_neighborsAndOrientations', '-v7.3');
end
% STAGE 2 ---- DEMIX ONCE ONLY (NEURONS MAY SEEM TO OVERLAP AT THE RESOLUTION OF LIGHT MICROSCOPY)
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0; for kk = 1:size(bbVol, 4); rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack; end; clear rawStack;
voxelCount                                                        = prod(stackSize);
detcov                                                            = zeros(1,numel(svCells));
for kk=1:numel(svCells)
  if numel(svCells{kk})>1; tmp=[]; for dd=1:size(bbVol, 4); tmp=[tmp bbVol(svCells{kk} + (dd-1)*voxelCount)]; end; tmp=tmp./repmat(sqrt(sum(tmp.^2,2)),1,size(tmp,2)); detcov(kk)=det(cov(tmp)); end;
end
[square_sAff, svMeans, svCells, voxelCounts, ~, ~, dups]          = demixSupervoxels5(                  square_sAff, svMeans, svCells, voxelCounts, detcov, opts.demix, 1, 1);
dupCount                                                          = size(dups,1);
thisFN                                                            = [opts.saveFileName 'stage2.mat'];
save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'dupCount', 'counts_neighborsAndOrientations', '-v7.3');
clear bbVol; % BY SPARSIFYING BBVOL, A VERY MEMORY EFFICIENT IMPLEMENTATION SHOULD BE ACHIEVED
% STAGE 3 ---- MERGE NEIGHBORING SUPERVOXELS WITH SIMILAR COLORS AND USING THE HEURISTICS OF STAGE 1
overClusteringFactor                                              = opts.kmeansMerging.overClusteringFactor;
while true
  colorTriplets                                                   = nchoosek(1:size(svMeans, 2), 3);
  allLUV                                                          = zeros(size(svMeans, 1), size(colorTriplets, 1)*3);
  for kk = 1:size(colorTriplets, 1)
    allLUV(:, (kk-1)*3+1:kk*3)                                    = rgb2luv(svMeans(:, colorTriplets(kk, :))')';
  end
  [coeff,score,latent]                                            = pca(allLUV);
  colorsForClustering                                             = score(:, 1:size(svMeans, 2));
  index                                                           = kmeans(colorsForClustering, round(opts.kmeansMerging.clusterCount*overClusteringFactor), 'MaxIter', 1000, 'Replicates', 100);
  binsaff_18n                                                     = (square_sAff>1/(sqrt(2)+1e-4));
  C                                                               = zeros(size(index));
  for kk = 1:max(index)
    theseSVs                                                      = find(index==kk);
    thisLUV                                                       = colorsForClustering(theseSVs, :);
    bincaff                                                       = (squareform(pdist(thisLUV))<opts.kmeansMerging.maxColorDistance);
    if isempty(bincaff)
      bincaff                                                     = true;
    end
    [thisS, thisC]                                                = graphconncomp(binsaff_18n(theseSVs, theseSVs) & bincaff, 'Weak', true);
    C(theseSVs)                                                   = max(C) + thisC;
  end
  C(index==0)                                                     = max(C)+1:max(C)+nnz(index==0);
  [square_sAff, svMeans, svCells, voxelCounts]                    = mergeClusters(C, [], square_sAff, svMeans, svCells, voxelCounts);
  opts.mergeSingleNeighborSuperVoxels.maxSizeForSingleNeighborSVs = quantile(voxelCounts, 0.5);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = newMergeWRTneighborsAndOrientations(square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeWRTnAo.normFlag, stackSize, opts.mergeWRTnAo);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeSmallSuperVoxels(              square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeSmallSuperVoxels);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeSingleNeighborSuperVoxels(     square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeSingleNeighborSuperVoxels);
  [square_sAff, svMeans, svCells, voxelCounts, ~]                 = mergeCloseNeighborhoods(            square_sAff, svMeans, svCells, voxelCounts, 1, opts.mergeCloseNeighborhoods);
  newCount = numel(svCells); if newCount==oldCount; break; end;
  if newCount<oldCount+10; overClusteringFactor=0.9*overClusteringFactor; end;
  if overClusteringFactor<2; break; end;
  oldCount = newCount; counts_neighborsAndOrientations(end+1) = oldCount; disp(counts_neighborsAndOrientations);
  thisFN = [opts.saveFileName 'stage3_iter' num2str(numel(counts_neighborsAndOrientations)-1) '.mat'];
  save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'dupCount', 'counts_neighborsAndOrientations', '-v7.3');
end
% STAGE 4 - SPATIAL AFFINITY CALCULATION
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
superVoxelOpts.spatialDistanceCalculationOpts.upperBound = 15;
sAff                                                     = calculate_sAff(cc, boundaryVoxelsSub, superVoxelOpts.spatialDistanceCalculationOpts);
[ii_sAff, ~, ss_sAff]                                    = find(sAff);
yy                                                       = ceil( (2*cc-1 - sqrt((2*cc-1)^2-8*ii_sAff))/2 );
xx                                                       = ii_sAff - cc*yy + cc + yy.*(yy+1)/2;
square_sAff                                              = sparse(xx, yy, ss_sAff, cc, cc);
square_sAff                                              = square_sAff + transpose(square_sAff);
clear ii_sAff; clear ss_sAff; clear xx; clear yy; clear sAff;
thisFN = [opts.saveFileName '_sAff15.mat'];
save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'dupCount', 'counts_neighborsAndOrientations', '-v7.3');
