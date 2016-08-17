
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPERVOXEL GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
superVoxelOpts.filePreamble                                              = '~/bb/data/dawen/supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.005_isplit0.5_30_augmented0.1';
superVoxelOpts.brightnessThreshold                                       = 0.1;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound                 = 2;
superVoxelOpts.splitInconsistentSVopts.maxPerimeter                      = 0.5;
superVoxelOpts.splitInconsistentSVopts.connectivity                      = 26;
superVoxelOpts.splitInconsistentSVopts.subdivisionSizeThreshold          = 20;
superVoxelOpts.dataset                                                   = '~/bb/data/dcai_PVhippo1_bm4d_sigma500.mat';
superVoxelOpts.HMINTH26                                                  = 0.005;
supervoxelize3(superVoxelOpts);
supervoxelize2(superVoxelOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSERVATIVE MERGING OF SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergeOpts.loadFilename                                                   = '/home/uygar/bb/data/dawen/supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.005_isplit0.5_30_augmented0.1_aff.mat';
mergeOpts.saveFileName                                                   = '/home/uygar/bb/data/dawen/merged_supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.005_isplit0.5_30_augmented0.1_minmax_demixFirst40_smallLast_maxCdist10';
mergeOpts.zAnisotropy                                                    = 3;
mergeOpts.demix.maxSimilarNeighborNormLUVDist                            = 40; % * sqrt(size(svMeans, 2)/4);
mergeOpts.demix.minImprovementFactor                                     = 3;
mergeOpts.demix.maxSizeForDemixing                                       = 500;
mergeOpts.mergeSmallSuperVoxels.luvColorDistanceUpperBound               = 20;
mergeOpts.mergeSmallSuperVoxels.disconnectedSVsizeTh                     = 20;
mergeOpts.mergeSmallSuperVoxels.maxVoxColorDist                          = 0.5;
mergeOpts.mergeWRTnAo.sDist                                              = sqrt(3);
mergeOpts.mergeWRTnAo.minDotProduct                                      = 0.9659; % pi/12 % sqrt(3)/2;
mergeOpts.mergeWRTnAo.maxColorDist                                       = 10; % * sqrt(size(svMeans, 2)/4);
mergeOpts.mergeWRTnAo.normFlag                                           = true;
mergeOpts.mergeWRTnAo.maxVoxColorDist                                    = 0.5;
mergeOpts.mergeSingleNeighborSuperVoxels.maxVoxColorDist                 = 0.5;
mergeOpts.mergeCloseNeighborhoods.maxDistNormLUV                         = 10; % * sqrt(size(svMeans, 2)/4);
mergeOpts.mergeCloseNeighborhoods.maxVoxColorDist                        = 0.5;
mergeSupervoxels2(mergeOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION OF MERGED SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergedSvFileName                    = '/home/uygar/bb/data/dawen/merged_supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.005_isplit0.5_30_augmented0.1_minmax_demixFirst40_smallLast_maxCdist10_sAff10.mat';
graphData.c                         = 1e-3;
graphData.s                         = 0;
graphData.colorRadiusForPure        = 20;
graphData.minSizeForPure            = 30;
graphData.maxPerim                  = 0.4;
graphData.spatialNhoodRadius        = sqrt(11)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.maxAssumedSdist           = 3;
graphData.opts_irbleigs.K           = 39;
[index, graphData]                  = segmentImage2(mergedSvFileName, graphData);
tempQscript; disp(max(hardClusterStats.ari)); %writeProjectedSegmentationScript;



load /home/uygar/bb/data/dawen/merged_supervoxels_dcai_PVhippo1_bm4d_sigma500_ws0.005_isplit0.5_30_augmented0.1_minmax_sAff10.mat
superVoxelCells=svCells;
perims                                                    = zeros(1, numel(superVoxelCells));
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0;
for kk = 1:size(bbVol, 4)
  rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack;
end
clear rawStack;
normalizer                                                = sqrt(sum(bbVol.^2, 4));
for dd = 1:size(bbVol, 4)
  bbVol(:,:,:,dd) = bbVol(:,:,:,dd) ./ normalizer;
end; clear normalizer;
voxelCount=prod(stackSize);
for kk = 1:numel(superVoxelCells)
  tmp                                                     = zeros(numel(superVoxelCells{kk}), size(bbVol, 4));
  for dd = 1:size(bbVol, 4)
    tmp(:, dd)                                            = bbVol(superVoxelCells{kk} + (dd-1)*voxelCount);
  end
  for dd = 1:size(bbVol, 4)
    perims(kk)                                            = max(perims(kk), max(tmp(:,dd))-min(tmp(:,dd)));
  end
end
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0;
for kk = 1:size(bbVol, 4)
  rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack;
end
clear rawStack;
[sorted,idx]=sort(perims,'descend');
cd ~/bb
im=showIndividualSuperVoxels(bbVol,svCells,idx(1:30));figure;imshow(im(:,:,1:3),[])

