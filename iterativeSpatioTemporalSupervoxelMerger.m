function [svMeans, svCells, square_sAff] = iterativeSpatioTemporalSupervoxelMerger(svMeans, svCells, square_sAff, stackSize, opts_mCCWIC, opts_fkmeans)

cc                                             = size(svMeans, 1);
ccOLD                                          = 2*cc;
voxelCounts                                    = zeros(cc, 1);
for kk = 1:cc
  voxelCounts(kk)                              = numel(svCells{kk});
end
nhoods                                         = [6 18 26];
nhood                                          = 6;
nhoodPointer                                   = 1;
flag                                           = false;
ccLog                                          = cc;
osFactorLog                                    = opts_mCCWIC.overSegmentingFactor;

while opts_mCCWIC.overSegmentingFactor>opts_mCCWIC.minOverSegmentingFactor & ccOLD~=cc
  svMeansNorm                                  = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
  if size(svMeans,2)==3
    svMeansLUV                                 = rgb2luv(svMeans')';
    svMeansNormLUV                             = rgb2luv(svMeansNorm')';
  else
    svMeansLUV(:,1:3)                          = rgb2luv(svMeans(:,1:3)')';
    svMeansNormLUV(:,1:3)                      = rgb2luv(svMeansNorm(:,1:3)')';
    svMeansLUV(:,4)                            = svMeans(:,4)*50;
    svMeansNormLUV(:,4)                        = svMeansNorm(:,4)*50;
  end
  index                                        = snrAwareKmeans(round(opts_mCCWIC.overSegmentingFactor*opts_mCCWIC.K), svMeansNormLUV, opts_fkmeans);
  [square_sAff, svMeans, svCells, voxelCounts] = mergeConnectedComponentsWithinIndividualClusters(index, square_sAff, svMeans, svCells, voxelCounts, opts_mCCWIC, nhood);
  opts_fkmeans.weight                          = voxelCounts;
  ccOLD                                        = cc;
  cc                                           = numel(svCells);
  ccLog(end+1)                                 = cc;
  osFactorLog(end+1)                           = opts_mCCWIC.overSegmentingFactor;
  if ccOLD/cc < opts_mCCWIC.fastRatio
    if flag
      opts_mCCWIC.overSegmentingFactor         = opts_mCCWIC.overSegmentingFactor * opts_mCCWIC.decayFactor;
      flag                                     = false;
    end
    nhoodPointer                               = rem(nhoodPointer, numel(nhoodPointer)) + 1;
    nhood                                      = nhood(nhoodPointer);
    if nhoodPointer == numel(nhoodPointer)
      flag                                     = true;
    end
  end
  disp(cc)
  if rem(numel(ccLog), 10) == 0
    save -v7.3 data/iSTSM_sv_dcai2_sigma2000_hmin26_0.008_26n_split5e-11_sCR_aff.mat ccLog osFactorLog index opts_fkmeans opts_mCCWIC square_sAff stackSize superVoxelOpts svCells svMeans voxelCounts
  end
end
