function [index, scoreInf, score1] = calculateNormalizedCuts(gD)

cc                                  = size(gD.colorData, 1);
opts_fkmeans.weight                 = sqrt(gD.svSizes);
opts_fkmeans.careful                = true;
opts_fkmeans.maxIter                = 1000;
opts_irbleigs                       = gD.opts_irbleigs;
opts_recon.cFactor1                 = gD.c1 * sqrt(3/size(gD.colorData, 2));
opts_recon.cFactor2                 = gD.c2 * sqrt(3/size(gD.colorData, 2));
opts_recon.color2SizeFactor         = gD.c2sf * median(gD.svSizes);
opts_recon.avgColorNeighborCount    = round(gD.nC * cc/opts_irbleigs.K);
opts_recon.sFactor1                 = gD.s1;
opts_recon.sFactor2                 = gD.s2;
opts_recon.sigma                    = gD.sig;

index                               = zeros(cc, 1);
[affinityMatrix, weightVector]      = calculateGraphAffinity(opts_recon, gD);
DD                                  = sparse(1:cc, 1:cc, 1./sqrt(sum(affinityMatrix)+eps));
[topFewEigenvectors, ss, PRGINF]    = irbleigs(opts_recon.sigma*DD + DD*affinityMatrix*DD, opts_irbleigs);
for ii=1:size(topFewEigenvectors,1)
  topFewEigenvectors(ii,:)          = topFewEigenvectors(ii,:) / norm(topFewEigenvectors(ii,:));
end
distortion                          = 1e60;
%allCentroids                       = zeros(opts_irbleigs.K, size(gD.colorData,2));
for kk = 1:50
  if isfield(gD,'ccIDsOfSVs')
    [myIndex, myCentroids, myDisto] = ss_fkmeans(topFewEigenvectors, opts_irbleigs.K, opts_fkmeans, gD.ccIDsOfSVs);
  else
    [myIndex, myCentroids, myDisto] = simple_fkmeans(topFewEigenvectors, opts_irbleigs.K, opts_fkmeans);
  end
  if sum(myDisto)<distortion
    distortion                      = sum(myDisto);
    index=myIndex;
%    allCentroids=myCentroids;
  end
end

clusterCount                        = max(index);
detcov                              = zeros(1, clusterCount);
for kk = 1:clusterCount
  tmp                               = gD.colorData(index==kk, :);
  tmp                               = tmp ./ repmat(sqrt(sum(tmp.^2,2)),1,size(tmp,2));
  detcov(kk)                        = det(cov(tmp));
end
scoreInf                            = max(detcov);
score1                              = sum(detcov);
