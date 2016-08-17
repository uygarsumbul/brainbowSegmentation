
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPERVOXEL GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
superVoxelOpts.filePreamble                                              = '~/bb/data/dawen/supervoxels_dcai_brainbow4channel_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1'; % PVhippo1_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1';
superVoxelOpts.brightnessThreshold                                       = 0.1;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound                 = 2;
superVoxelOpts.splitHardSVopts.detThreshold                              = 5e-11;
superVoxelOpts.splitHardSVopts.connectivity                              = 26;
superVoxelOpts.splitHardSVopts.subdivisionSizeThreshold                  = 30;
superVoxelOpts.dataset                                                   = '~/bb/data/dcai_brainbow4channel_bm4d_sigma500.mat'; % dcai_PVhippo1_bm4d_sigma500.mat';
superVoxelOpts.HMINTH26                                                  = 0.005;
supervoxelize(superVoxelOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSERVATIVE MERGING OF SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergeOpts.loadFilename                                                   = '/home/uygar/bb/data/dawen/supervoxels_dcai_brainbow4channel_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1_aff.mat'; % PVhippo1_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1_aff.mat';
mergeOpts.saveFileName                                                   = '/home/uygar/bb/data/dawen/merged_supervoxels_dcai_brainbow4channel_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1'; % PVhippo1_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1';
mergeOpts.zAnisotropy                                                    = 3;
mergeOpts.demix.maxSimilarNeighborNormLUVDist                            = 30; % * sqrt(size(svMeans, 2)/4);
mergeOpts.demix.minImprovementFactor                                     = 3;
mergeOpts.demix.maxSizeForDemixing                                       = 500;
mergeOpts.mergeSmallSuperVoxels.luvColorDistanceUpperBound               = 30;
mergeOpts.mergeSmallSuperVoxels.disconnectedSVsizeTh                     = 20;
mergeOpts.mergeWRTnAo.sDist                                              = sqrt(3);
mergeOpts.mergeWRTnAo.minDotProduct                                      = 1/sqrt(2);
mergeOpts.mergeWRTnAo.maxColorDist                                       = 6; % * sqrt(size(svMeans, 2)/4);
mergeOpts.mergeWRTnAo.normFlag                                           = true;
mergeOpts.mergeCloseNeighborhoods.maxDistNormLUV                         = 6; % * sqrt(size(svMeans, 2)/4);
mergeOpts.kmeansMerging.clusterCount                                     = 40;
mergeOpts.kmeansMerging.overClusteringFactor                             = 5;
mergeOpts.kmeansMerging.maxColorDistance                                 = 6; % * sqrt(size(svMeans, 2)/4);
mergeSupervoxels(mergeOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION OF MERGED SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergedSvFileName                    = '/home/uygar/bb/data/dawen/merged_supervoxels_dcai_PVhippo1_bm4d_sigma500_ws0.005_split5e-11_30_augmented0.1_sAff15.mat';
graphData.c                         = 1e-3;
graphData.s                         = 4.33e-2;
graphData.colorRadiusForPure        = 20;
graphData.minSizeForPure            = 50;
graphData.maxDetCovForPure          = 1e-10;
graphData.spatialNhoodRadius        = sqrt(11)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.maxAssumedSdist           = 1;
graphData.opts_irbleigs.K           = 39;
index                               = segmentImage(mergedSvFileName, graphData);
writeProjectedSegmentationScript;
