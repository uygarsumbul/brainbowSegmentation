function [svCells, origIndex] = scMergerOnSubgraph(osFactor, minOsFactor, K_vdpgm, opts, square_sAff, svMeans, svCells, voxelCounts)

disp([K_vdpgm*10 numel(svCells)])

allConnTh                                                    = [0.9 0.7 0.51]; % corresponding neighborhoods: 6, 18, 26                                                                             
connThPointer                                                = 1;
connTh                                                       = allConnTh(connThPointer);
cc                                                           = size(svMeans, 1);
ccOLD                                                        = cc + 1;
normFlag                                                     = true;
origIndex                                                    = 1:cc;
while osFactor>minOsFactor & ccOLD~=cc
  if normFlag
    colorData                                                = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
    if size(svMeans,2)==3
      colorData                                              = rgb2luv(colorData')';
    else
      colorData                                              = [rgb2luv(colorData(:,1:3)')' colorData(:,4)*50];
    end
  else
    if size(svMeans,2)==3
      colorData                                              = rgb2luv(svMeans')';
    else
      colorData                                              = [rgb2luv(svMeans(:,1:3)')' svMeans(:,4)*50];
    end
  end
  doSpectralSplit                                            = (osFactor<=opts.maxOsFactorForSC);
  opts.fkmeans.weight                                        = voxelCounts;
  index                                                      = snrAwareKmeans(round(osFactor*K_vdpgm), colorData, opts.fkmeans);
  [square_sAff, svMeans, svCells, voxelCounts, idxTransform] = mergeConnComp3(index, square_sAff, svMeans, svCells, voxelCounts, opts, connTh, doSpectralSplit);
  origIndex                                                  = idxTransform(origIndex);
  ccOLD                                                      = cc;
  cc                                                         = numel(svCells);
  if ccOLD/cc < opts.fastRatio
    if connThPointer == numel(allConnTh)
      osFactor                                               = max(minOsFactor-eps, osFactor * opts.decayFactor);
    end
    connThPointer                                            = rem(connThPointer, numel(allConnTh)) + 1;
    connTh                                                   = allConnTh(connThPointer);
  end
  normFlag                                                   = ~normFlag;
end
