function [index, scoreInf, score1] = calculateNormalizedCuts_mem(graphData)

cc                                  = size(graphData.colorData, 1);
opts_fkmeans.weight                 = sqrt(graphData.svSizes);
opts_fkmeans.careful                = true;
opts_fkmeans.maxIter                = 1000;
opts_irbleigs                       = graphData.opts_irbleigs;
opts_recon.cFactor1                 = graphData.c1 * sqrt(3/size(graphData.colorData, 2));
opts_recon.cFactor2                 = graphData.c2 * sqrt(3/size(graphData.colorData, 2));
opts_recon.color2SizeFactor         = graphData.c2sf * mean(graphData.svSizes);
opts_recon.avgColorNeighborCount    = round(graphData.nC * cc/opts_irbleigs.K);
opts_recon.sFactor1                 = graphData.s1;
opts_recon.sFactor2                 = graphData.s2;
opts_recon.sigma                    = graphData.sig;
index                               = zeros(cc, 1);

if ~isfield(graphData,'seedSets')   || isempty(graphData.seedSets)  ; seedSets   = cell(0);                                else; seedSets   = graphData.seedSets;   end;
if ~isfield(graphData,'ccIDsOfSVs') || isempty(graphData.ccIDsOfSVs); ccIDsOfSVs = zeros(size(graphData.colorData, 1), 1); else; ccIDsOfSVs = graphData.ccIDsOfSVs; end;
colorData                                    = graphData.colorData;
stackSize                                    = graphData.stackSize;
graphData.svSizes                            = transpose(graphData.svSizes(:));
localTh                                      = 30;
minNeighborCount                             = round(0.0005 * numel(graphData.svSizes)/graphData.opts_irbleigs.K); %5;
%% CALCULATE SQUAREFORM AFFINITY MATRIX AND RESTRICT VARIABLES TO THE ACTIVE NODE SUBSET
voxelCount                                   = prod(stackSize);
xx_cAff                                      = zeros(round((graphData.nC/opts_irbleigs.K)*cc*(cc-1)/2), 1);
yy_cAff                                      = zeros(round((graphData.nC/opts_irbleigs.K)*cc*(cc-1)/2), 1);
ss_cAff                                      = zeros(round((graphData.nC/opts_irbleigs.K)*cc*(cc-1)/2), 1);
cIdx                                         = 1;
Mdl                                          = KDTreeSearcher(colorData);
luvColorDistanceUpperBound                   = selectColorDistance(Mdl, colorData, opts_recon);
disp(luvColorDistanceUpperBound)
%% GENERATE CANNOT-CONNECT SETS FROM SEED SETS
cantConnect                                  = cell(numel(seedSets), 1);
for kk1 = 1:numel(seedSets)
  for kk2 = 1:numel(seedSets)
    if kk2 ~= kk1
      cantConnect{kk1}                       = [cantConnect{kk1}; seedSets{kk2}];
    end
  end
end
%% LOOP OVER INDIVIDUAL SUPERVOXELS
for kk1 = 1:cc
  % COLOR NEIGHBORS
  [colorNeighborsCell, Dcell]                = rangesearch(Mdl, colorData(kk1, :), luvColorDistanceUpperBound);
  % MAKE SURE SOME COLOR NEIGHBORS ARE SELECTED
  tmpBound                                   = luvColorDistanceUpperBound;
  tmpBound                                   = tmpBound * 1.05;
  while numel(Dcell{1})<minNeighborCount+1
    [colorNeighborsCell, Dcell]              = rangesearch(Mdl, colorData(kk1, :), tmpBound);
    tmpBound                                 = tmpBound * 1.05;
  end
  % REMOVE SELF
  D                                          = Dcell{1}(2:end);
  colorNeighbors                             = colorNeighborsCell{1}(2:end);
  % CORRESPONDING SPATIAL AFFINITIES
  this_sAff                                  = graphData.square_sAff(kk1, colorNeighbors);
  % RETAIN NEIGHBORS THAT ARE NOT FAR AWAY IN SPACE (SIZE-WEIGHTED SPATIAL DISTANCE)
  harmmeanVoxCounts                          = harmmean([graphData.svSizes(kk1)*ones(1, numel(colorNeighbors)); graphData.svSizes(colorNeighbors)]);
  tmp                                        = find(this_sAff+1/localTh >= min(0.5, opts_recon.color2SizeFactor ./ harmmeanVoxCounts));
  if numel(tmp)<minNeighborCount
    tmp2                                     = setdiff(1:numel(D), tmp, 'stable');
    tmp                                      = [tmp tmp2(1:max(0, minNeighborCount-numel(tmp)))];
  end
  D                                          = D(tmp);
  colorNeighbors                             = colorNeighbors(tmp);
  D                                          = exp( (-opts_recon.cFactor1 - opts_recon.cFactor2*log(harmmeanVoxCounts(tmp))) .* (D.^2) );
  D                                          = D .* (opts_recon.sFactor1 + (1-opts_recon.sFactor1)*exp( -opts_recon.sFactor2*(-1+1./this_sAff(tmp)).^2 ));
  if ccIDsOfSVs(kk1)>0
    % CANNOT-CONNECT
    toRemove                                 = ismember(colorNeighbors, cantConnect{ccIDsOfSVs(kk1)});
    colorNeighbors(toRemove)                 = [];
    D(toRemove)                              = [];
    % MUST-CONNECT
    thisSeed                                 = setdiff(seedSets{ccIDsOfSVs(kk1)}, kk1, 'stable');
    [~, ia, ib]                              = intersect(colorNeighbors, thisSeed);
    D(ia)                                    = 1;
    remaining                                = setdiff(1:numel(thisSeed), ib);
    colorNeighbors                           = [colorNeighbors thisSeed(remaining)'];
    D                                        = [D ones(1, numel(remaining))];
  end
  xx_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = kk1;
  yy_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = colorNeighbors;
  ss_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = D;
  cIdx                                       = cIdx+numel(colorNeighbors);
end
unused                                       = find(~xx_cAff, 1, 'first');
if ~isempty(unused)
  ss_cAff(unused:end)                        = [];
  xx_cAff(unused:end)                        = [];
  yy_cAff(unused:end)                        = [];
end
affinityMatrix                               = sparse(xx_cAff, yy_cAff, ss_cAff, cc, cc);
affinityMatrix                               = max(affinityMatrix, transpose(affinityMatrix));
weightVector                                 = full(sum(affinityMatrix, 2));


DD                                  = sparse(1:cc, 1:cc, 1./sqrt(sum(affinityMatrix)+eps));
[topFewEigenvectors, ss, PRGINF]    = irbleigs(opts_recon.sigma*DD + DD*affinityMatrix*DD, opts_irbleigs);
for ii=1:size(topFewEigenvectors,1)
  topFewEigenvectors(ii,:)          = topFewEigenvectors(ii,:) / norm(topFewEigenvectors(ii,:));
end
distortion                          = 1e60;
%allCentroids                       = zeros(opts_irbleigs.K, size(graphData.colorData,2));
for kk = 1:50
  if isfield(graphData,'ccIDsOfSVs')
    [myIndex, myCentroids, myDisto] = ss_fkmeans(topFewEigenvectors, opts_irbleigs.K, opts_fkmeans, graphData.ccIDsOfSVs);
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
  tmp                               = graphData.colorData(index==kk, :);
  tmp                               = tmp ./ repmat(sqrt(sum(tmp.^2,2)),1,size(tmp,2));
  detcov(kk)                        = det(cov(tmp));
end
scoreInf                            = max(detcov);
score1                              = sum(detcov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                        
function luvColorDistanceUpperBound = selectColorDistance(Mdl, colorData, opts_recon)
cc                             = size(colorData, 1);
samples                        = randi(cc, round(0.01*cc), 1);
lower                          = 0;
upper                          = 80;
lowerCount                     = 0;
upperCount                     = 0;
countEpsilon                   = 1;
avgColorNeighborCount          = opts_recon.avgColorNeighborCount;
for kk = 1:numel(samples)
  [tmpCell, Dcell]             = rangesearch(Mdl, colorData(samples(kk), :), lower);
  lowerCount                   = lowerCount + numel(tmpCell{1}); % nnz(exp( -opts_recon.cFactor * Dcell{1} )>=opts_recon.sparsifyingFactor );
  [tmpCell, Dcell]             = rangesearch(Mdl, colorData(samples(kk), :), upper);
  upperCount                   = upperCount + numel(tmpCell{1}); % nnz(exp( -opts_recon.cFactor * Dcell{1} )>=opts_recon.sparsifyingFactor );
end
lowerCount                     = lowerCount/numel(samples) - 1;
upperCount                     = upperCount/numel(samples) - 1;
if abs(lowerCount-avgColorNeighborCount)<countEpsilon
  luvColorDistanceUpperBound   = lower;
  return;
end
if abs(upperCount-avgColorNeighborCount)<countEpsilon
  luvColorDistanceUpperBound   = upper;
  return;
end
distanceList                   = [lower upper];
countList                      = [lowerCount upperCount];
while true
  pos                          = find(countList>avgColorNeighborCount, 1);
  if pos>1
    newDist                    = mean(distanceList(pos-1:pos));
  else
    newDist                    = distanceList(1)/2;
  end
  thisCount                    = 0;
  for kk = 1:numel(samples)
    [tmpCell, Dcell]           = rangesearch(Mdl, colorData(samples(kk), :), newDist);
    thisCount                  = thisCount + numel(tmpCell{1}); % nnz(exp( -opts_recon.cFactor * Dcell{1} )>=opts_recon.sparsifyingFactor );
  end
  thisCount                    = thisCount/numel(samples) - 1;
  if abs(thisCount-avgColorNeighborCount)<countEpsilon
    luvColorDistanceUpperBound = newDist;
    break;
  else
    if pos>1
      distanceList             = [distanceList(1:pos-1) newDist distanceList(pos:end)];
      countList                = [countList(1:pos-1) thisCount countList(pos:end)];
    else
      distanceList             = [newDist distanceList];
      countList                = [thisCount countList];
    end
  end
end
