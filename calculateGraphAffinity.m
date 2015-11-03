function [affinityMatrix, weightVector] = calculateGraphAffinity(opts_recon, graphData)

if ~isfield(graphData,'seedSets')   || isempty(graphData.seedSets)  ; seedSets   = cell(0);                                else; seedSets   = graphData.seedSets;   end;
if ~isfield(graphData,'ccIDsOfSVs') || isempty(graphData.ccIDsOfSVs); ccIDsOfSVs = zeros(size(graphData.colorData, 1), 1); else; ccIDsOfSVs = graphData.ccIDsOfSVs; end;
if ~isfield(graphData,'svSubset')   || isempty(graphData.svSubset)  ; svSubset   = 1:size(graphData.colorData, 1);         else; svSubset   = graphData.svSubset;   end;
colorData                                    = graphData.colorData;
stackSize                                    = graphData.stackSize;
svSizes                                      = graphData.svSizes;
square_sAff                                  = graphData.square_sAff;
svSubset                                     = sort(svSubset);
mapping                                      = 1:numel(svSubset);
localTh                                      = 30;
minNeighborCount                             = 5;
%% CALCULATE SQUAREFORM AFFINITY MATRIX AND RESTRICT VARIABLES TO THE ACTIVE NODE SUBSET
ccIDsOfSVs                                   = ccIDsOfSVs(svSubset);
cc                                           = size(colorData, 1); % reassign below
voxelCount                                   = prod(stackSize);
square_sAff                                  = square_sAff(svSubset, svSubset);
cc                                           = numel(svSubset);
svSizes                                      = svSizes(svSubset);
svSizes                                      = transpose(svSizes(:));
ccIDsOfSVs                                   = ccIDsOfSVs(svSubset);
for kk = 1:numel(seedSets)
  seedSets{kk}                               = mapping(intersect(seedSets{kk}, svSubset, 'stable'));
end

xx_cAff                                      = zeros(round(0.05*cc*(cc-1)/2), 1);
yy_cAff                                      = zeros(round(0.05*cc*(cc-1)/2), 1);
ss_cAff                                      = zeros(round(0.05*cc*(cc-1)/2), 1);
cIdx                                         = 1;
Mdl                                          = KDTreeSearcher(colorData);
luvColorDistanceUpperBound                   = selectColorDistance(Mdl, colorData, opts_recon);
disp(luvColorDistanceUpperBound)
%% GENERATE CANNOT-CONNECT SETS FROM SEED SETS
cantConnect                                  = cell(numel(seedSets), 1);
for kk1 = 1:numel(seedSets)
  for kk2 = 1:numel(seedSets)
    if kk2 ~= kk1
      cantConnect{kk1}                       = [cantConnect{kk1} seedSets{kk2}];
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
  this_sAff                                  = square_sAff(kk1, colorNeighbors);
  % RETAIN NEIGHBORS THAT ARE NOT FAR AWAY IN SPACE (SIZE-WEIGHTED SPATIAL DISTANCE)
  harmmeanVoxCounts                          = harmmean([svSizes(kk1)*ones(1, numel(colorNeighbors)); svSizes(colorNeighbors)]);
  tmp                                        = find(this_sAff+1/localTh >= min(0.5, opts_recon.color2SizeFactor ./ harmmeanVoxCounts));
  % ENFORCE CONNECTIVITY FOR NEIGHBORS WITHIN THE SPATIAL NEIGHBORHOOD
%  tmpLocal                                   = intersect(tmp, find(this_sAff >= 1/localTh), 'stable');
%  tmpFar                                     = setdiff(tmp, tmpLocal);
%  localGraph                                 = square_sAff([colorNeighbors(tmpLocal) kk1], [colorNeighbors(tmpLocal) kk1]);
%  localGraph                                 = localGraph .* (localGraph>1/(sqrt(1)+1e-6));    % 6-CONNECTIVITY!!!!!!!!!!!!!!
%  [S, C]                                     = graphconncomp(localGraph);
%  tmpLocal                                   = tmpLocal(find(C(1:end-1)==C(end)));
%  tmp                                        = [tmpLocal tmpFar];
  % ENSURE THAT AT LEAST TWO COLOR NEIGHBORS ARE CHOSEN
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
    colorNeighbors                           = [colorNeighbors thisSeed(remaining)];
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
