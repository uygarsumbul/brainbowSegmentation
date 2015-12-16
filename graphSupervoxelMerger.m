function [svMeans, svCells, square_sAff, origIndex] = graphSupervoxelMerger(svMeans, svCells, square_sAff, opts)

% opts. -- fkmeans, superVoxelOpts, stackSize, maxOverSegmentingFactor, minOverSegmentingFactor, maxOsFactorForSC, decayFactor, maxSubgraphSize, maxColorDifferenceRatio, disconnectedSVsizeTh
cc                                             = size(svMeans, 1);
ccOLD                                          = 2*cc;
osFactor                                       = opts.maxOverSegmentingFactor;
voxelCounts                                    = cellfun(@numel, svCells)';
ccLog                                          = cc;
osFactorLog                                    = osFactor;
origIndex                                      = 1:cc;
while osFactor>opts.minOverSegmentingFactor % & ccOLD~=cc
  subgraphNodes                                = randomlyPartitionBinaryGraph(square_sAff > 0.9, opts.maxSubgraphSize);
  subgraph_svCells                             = cell(1, numel(subgraphNodes));
  subgraph_indices                             = cell(1, numel(subgraphNodes));
  nextOsFactor                                 = opts.outerDecayFactor*osFactor;
  if size(svMeans,2)==3
    svMeansLUV                                 = rgb2luv(svMeans')';
  else
    svMeansLUV                                 = [rgb2luv(svMeans(:,1:3)')' svMeans(:,4)*50];
  end
  parfor sg = 1:numel(subgraphNodes)
    nn                                         = subgraphNodes{sg};
    K_vdpgm                                    = vdpgm_countEstimate(svMeansLUV(nn, :));
    [sg_svCells, sg_svIndex]                   = scMergerOnSubgraph(osFactor, nextOsFactor, K_vdpgm, opts, square_sAff(nn, nn), svMeans(nn, :), svCells(nn), voxelCounts(nn));
    subgraph_svCells{sg}                       = sg_svCells;
    subgraph_indices{sg}                       = sg_svIndex;
  end
  [square_sAff, svMeans, svCells, voxelCounts, idxTransform] = mergeSubgraphs2(square_sAff, svMeans, voxelCounts, opts.disconnectedSVsizeTh, subgraph_svCells, subgraph_indices, subgraphNodes);
  origIndex                                                  = idxTransform(origIndex);
  ccOLD                                        = cc;
  cc                                           = numel(svCells);
  ccLog(end+1)                                 = cc;
  osFactorLog(end+1)                           = osFactor;
  disp([cc osFactorLog(end)])
  osFactor                                     = nextOsFactor;
  if isfield(opts, 'saveFileName')
    fileName = [opts.saveFileName '_iter' num2str(numel(ccLog)-1) '.mat'];
  else
    fileName = ['/home/uygar/bb/data/iSTSM_sv_dcai2_sigma2000_hmin26_0.008_26n_split5e-11_sCR_iter' num2str(numel(ccLog)-1) '.mat'];
  end
  save(fileName, 'ccLog', 'osFactorLog', 'opts', 'square_sAff', 'svCells', 'svMeans', 'voxelCounts', 'origIndex', '-v7.3');
end

