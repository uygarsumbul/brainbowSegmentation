function [index, graphData] = segmentImage(fileName, graphData)

load(fileName);
voxelCount                          = prod(stackSize);
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0; for kk = 1:size(bbVol, 4); rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack; end; clear rawStack;
graphData.svCells                   = svCells;
graphData.sig                       = 0;
opts_fkmeans.careful                = true;
opts_fkmeans.maxIter                = 20000;
opts_recon.sigma                    = graphData.sig;
svMeansNorm                         = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
allTriplets                         = nchoosek(1:size(svMeans, 2), 3);
allColors                           = zeros(size(svMeans, 1), 3*size(allTriplets, 1));
allColorsNorm                       = allColors;
svMeansNorm                         = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
for kk = 1:size(allTriplets, 1)
  allColors(:, 3*kk-2:3*kk)         = rgb2luv(svMeans(:, allTriplets(kk, :))')';
  allColorsNorm(:, 3*kk-2:3*kk)     = rgb2luv(svMeansNorm(:, allTriplets(kk, :))')';
end
[~, score, ~]                       = pca(allColors);
graphData.colorsForDistal           = score(:, 1:size(svMeans, 2));
[~, scoreNorm, ~]                   = pca(allColorsNorm);
graphData.colorsForProximal         = score(:, 1:size(svMeans, 2));
graphData.square_sAff               = square_sAff;
graphData.svSizes                   = voxelCounts;
graphData.smallSVs                  = 1:numel(svCells);

normalizer                          = sqrt(sum(bbVol.^2, 4));
perims                              = zeros(1,numel(svCells));
for dd = 1:size(bbVol, 4)
  bbVol(:,:,:,dd)                   = bbVol(:,:,:,dd) ./ normalizer;
end; clear normalizer;
for kk = 1:numel(svCells)
  tmp                               = zeros(numel(svCells{kk}), size(bbVol, 4));
  for dd = 1:size(bbVol, 4)
    tmp(:, dd)                      = bbVol(svCells{kk} + (dd-1)*voxelCount);
  end
  for dd = 1:size(bbVol, 4)
    perims(kk)                      = max(perims(kk), max(tmp(:,dd))-min(tmp(:,dd)));
  end
end
graphData.perims                    = perims;

affinityMatrix                      = generateAffinityMatrix2(graphData);
nodeDegrees                         = sum(affinityMatrix);
toRemove                            = find(nodeDegrees==0);
toKeep                              = setdiff(1:numel(svCells), toRemove); if ~isempty(toRemove); disp(['toRemove count: ' num2str(numel(toRemove))]); end;
thisAffinityMatrix                  = affinityMatrix;
thisAffinityMatrix(toRemove, :)     = [];
thisAffinityMatrix(:, toRemove)     = [];
opts_fkmeans.weight                 = graphData.svSizes(toKeep);
cc                                  = size(thisAffinityMatrix, 1);
weightVector                        = full(sum(thisAffinityMatrix, 2));
DD                                  = sparse(1:cc, 1:cc, 1./sqrt(sum(thisAffinityMatrix)+eps));
[topFewEigenvectors, ss, PRGINF]    = irbleigs(opts_recon.sigma*DD + DD*thisAffinityMatrix*DD, graphData.opts_irbleigs);
for ii=1:size(topFewEigenvectors,1)
  topFewEigenvectors(ii,:)          = topFewEigenvectors(ii,:) / norm(topFewEigenvectors(ii,:));
end
warning off;
distortion                          = 1e60;
for kk = 1:1000
  initialLabels                     = kmeans(graphData.colorsForDistal(toKeep, :), graphData.opts_irbleigs.K,'MaxIter',2);
  [myIndex, myCentroids, myDisto]   = colorInitializedWeightedKmeans(topFewEigenvectors, graphData.opts_irbleigs.K, initialLabels, opts_fkmeans);
  if sum(myDisto)<distortion;  distortion = sum(myDisto); index_minCut=myIndex; end;
end
index                               = zeros(1, numel(svCells)); index(toKeep) = index_minCut;
