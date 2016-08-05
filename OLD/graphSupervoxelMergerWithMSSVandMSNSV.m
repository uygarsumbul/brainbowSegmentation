function [svMeans, svCells, square_sAff, origIndex] = graphSupervoxelMergerWithMSSV(svMeans, origsvMeansLUV, origIndex, svCells, square_sAff, opts)

% opts. -- fkmeans, superVoxelOpts, stackSize, maxOverSegmentingFactor, minOverSegmentingFactor, maxOsFactorForSC, decayFactor, maxSubgraphSize, maxColorDifferenceRatio, disconnectedSVsizeTh
osFactor                                                     = opts.maxOverSegmentingFactor;
voxelCounts                                                  = cellfun(@numel, svCells)';
ccLog                                                        = size(svMeans, 1);
osFactorLog                                                  = osFactor;
origDisconnectedSVsizeTh                                     = opts.disconnectedSVsizeTh;
while osFactor>opts.minOverSegmentingFactor
  % MERGE BASED ON SPACE-COLOR PROXIMITY
  subgraphNodes                                              = randomlyPartitionBinaryGraph(square_sAff > 1/(sqrt(3)+1e-4), opts.maxSubgraphSize);
  subgraph_svCells                                           = cell(1, numel(subgraphNodes));
  subgraph_indices                                           = cell(1, numel(subgraphNodes));
  nextOsFactor                                               = opts.outerDecayFactor*osFactor;
  if size(svMeans,2)==3
    svMeansLUV                                               = rgb2luv(svMeans')';
  else
    svMeansLUV                                               = [rgb2luv(svMeans(:,1:3)')' svMeans(:,4)*50];
  end
  parfor sg = 1:numel(subgraphNodes)
    nn                                                       = subgraphNodes{sg};
    K_vdpgm                                                  = vdpgm_countEstimate(origsvMeansLUV(ismember(origIndex, nn), :));
    [sg_svCells, sg_svIndex]                                 = scMergerOnSubgraph(osFactor, nextOsFactor, K_vdpgm, opts, square_sAff(nn, nn), svMeans(nn, :), svCells(nn), voxelCounts(nn));
    subgraph_svCells{sg}                                     = sg_svCells;
    subgraph_indices{sg}                                     = sg_svIndex;
  end
  [square_sAff, svMeans, svCells, voxelCounts, idxTransform] = mergeSubgraphs2(square_sAff, svMeans, voxelCounts, opts.disconnectedSVsizeTh, subgraph_svCells, subgraph_indices, subgraphNodes);
  origIndex                                                  = idxTransform(origIndex);
  ccLog(end+1)                                               = numel(svCells);
  osFactorLog(end+1)                                         = osFactor;
  disp([ccLog(end) osFactorLog(end)])
  % MERGE SINGLE NEIGHBORS
  [square_sAff, svMeans, svCells, voxelCounts, origIndex] = mergeSingleNeighborSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, origIndex, opts);
  % MERGE BASED ON SIZE
  if osFactor<opts.maxOsFactorForSC
    opts.disconnectedSVsizeTh                                = round(mean(voxelCounts)/10);
    opts.luvColorDistanceUpperBound                          = 300/opts.disconnectedSVsizeTh;
    [square_sAff, svMeans, svCells, voxelCounts, origIndex]  = mergeSmallSuperVoxels(square_sAff, svMeans, svCells, voxelCounts, origIndex, opts);
    opts.disconnectedSVsizeTh                                = origDisconnectedSVsizeTh;
    ccLog(end+1)                                             = numel(svCells);
    osFactorLog(end+1)                                       = osFactor;
    disp([ccLog(end) osFactorLog(end)])
  end
  % UPDATE FOR NEXT ITERATION
  osFactor                                                   = nextOsFactor;
  opts.maxColorDistSizeProduct                               = min(opts.maxColorDistSizeProduct*opts.colorDistSizeProductMultiplier, opts.upperBoundForColorDistSizeProduct);
  fileName                                                   = [opts.saveFileName '_iter' num2str(numel(ccLog)-1) '.mat'];
  save(fileName, 'ccLog', 'osFactorLog', 'opts', 'square_sAff', 'svCells', 'svMeans', 'origsvMeansLUV', 'voxelCounts', 'origIndex', '-v7.3');
end

