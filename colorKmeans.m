function index = snrAwareKmeans(clusterCount, colorData, opts_fkmeans)

cc                                           = size(colorData, 1);
opts_irbleigs.K                              = clusterCount;

index                                        = zeros(cc, 1);
distortion                                   = 1e60;
allCentroids                                 = zeros(opts_irbleigs.K, size(colorData, 2));
for kk = 1:50
  [myIndex, myCentroids, myDistortion]       = fkmeans(colorData, opts_irbleigs.K, opts_fkmeans);
  if sum(myDistortion)<distortion
    distortion                               = sum(myDistortion);
    index=myIndex;
    allCentroids=myCentroids;
  end
end
