
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPERVOXEL GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
superVoxelOpts.dataset                                                   = '~/bb/data/dcai_brainbow4channel_bm4d_sigma500.mat'; % '~/bb/data/dcai_PVhippo1_bm4d_sigma500.mat';
superVoxelOpts.filePreamble                                              = '~/bb/data/dawen/supervoxels3_dcai_brainbow4channel_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1'; % supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.004_isplit0.4_20_augmented0.1';
superVoxelOpts.brightnessThreshold                                       = 0.1;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound                 = 2;
superVoxelOpts.splitInconsistentSVopts.maxPerimeter                      = 0.5;
superVoxelOpts.splitInconsistentSVopts.connectivity                      = 26;
superVoxelOpts.splitInconsistentSVopts.subdivisionSizeThreshold          = 20;
superVoxelOpts.HMINTH26                                                  = 0.006;
supervoxelize(superVoxelOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSERVATIVE MERGING OF SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergeOpts.loadFilename                                                   = '~/bb/data/dawen/supervoxels3_dcai_brainbow4channel_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1_aff.mat'; % PVhippo1_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1_aff.mat';
mergeOpts.saveFileName                                                   = '~/bb/data/dawen/merged_supervoxels3_dcai_brainbow4channel_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1_minmax_demixFirst50_5_smallLast_maxCdist10';
%PVhippo1_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1_minmax_demixFirst50_5_smallLast_maxCdist10';
mergeOpts.zAnisotropy                                                    = 3;
mergeOpts.demix.maxSimilarNeighborNormLUVDist                            = 50; % * sqrt(size(svMeans, 2)/4);
mergeOpts.demix.minImprovementFactor                                     = 5;
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
mergeSupervoxels(mergeOpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTATION OF MERGED SUPERVOXELS %%%%%%%%%%%%%%%%%%%%%%%%%
mergedSvFileName                    = '~/bb/data/dawen/merged_supervoxels3_dcai_brainbow4channel_bm4d_sigma500_ws0.006_isplit0.5_20_augmented0.1_minmax_demixFirst50_5_smallLast_maxCdist10_sAff10.mat';
%'~/bb/data/dawen/merged_supervoxels3_dcai_PVhippo1_bm4d_sigma500_ws0.004_isplit0.4_20_augmented0.1_minmax_demixFirst50_5_smallLast_maxCdist15_sAff10.mat';
graphData.c                         = 2e-3;
graphData.colorRadiusForPure        = 20;
graphData.minSizeForPure            = 100;
graphData.maxPerim                  = 0.4;
graphData.spatialNhoodRadius        = sqrt(9)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.opts_irbleigs.K           = 39;
[index, graphData]                  = segmentImage3(mergedSvFileName, graphData);
tempQscript; disp(max(hardClusterStats.ari)); %writeProjectedSegmentationScript;
