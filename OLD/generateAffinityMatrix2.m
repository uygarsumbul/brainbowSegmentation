function affinityMatrix = generateAffinityMatrix(graphData)

[cc, ch]                                     = size(graphData.colorsForDistal);
opts_irbleigs                                = graphData.opts_irbleigs;
opts_recon.cFactor                           = graphData.c * sqrt(3/ch);
opts_recon.sigma                             = graphData.sig;
index                                        = zeros(cc, 1);
if ~isfield(graphData,'seedSets')   || isempty(graphData.seedSets)  ; seedSets   = cell(0);      else; seedSets   = graphData.seedSets;   end;
if ~isfield(graphData,'ccIDsOfSVs') || isempty(graphData.ccIDsOfSVs); ccIDsOfSVs = zeros(cc, 1); else; ccIDsOfSVs = graphData.ccIDsOfSVs; end;
graphData.svSizes                            = transpose(graphData.svSizes(:));
%% CALCULATE SQUAREFORM AFFINITY MATRIX AND RESTRICT VARIABLES TO THE ACTIVE NODE SUBSET
xx_cAff                                      = zeros(round((1e-2/opts_irbleigs.K)*cc*(cc-1)/2), 1);
yy_cAff                                      = zeros(round((1e-2/opts_irbleigs.K)*cc*(cc-1)/2), 1);
ss_cAff                                      = zeros(round((1e-2/opts_irbleigs.K)*cc*(cc-1)/2), 1);
cIdx                                         = 1;
pure                                         = intersect(graphData.smallSVs, find(graphData.svSizes>graphData.minSizeForPure & graphData.perims<graphData.maxPerim)); % 30, 1e-10
Mdl                                          = KDTreeSearcher(graphData.colorsForDistal(pure, :));
%%%%% LUV-NORM COLOR NEIGHBORS WITHIN A RADIUS FOR SMALL SUPERVOXELS
for kk = 1:numel(graphData.smallSVs)
  kk1                                        = graphData.smallSVs(kk);
  colorNeighbors = [];
  D = [];
  if ismember(kk1, pure)
    [colorNeighbors, D]                      = rangesearch(Mdl,graphData.colorsForDistal(kk1, :),graphData.colorRadiusForPure*sqrt(ch/4));
    colorNeighbors                           = pure(colorNeighbors{1}(2:end));
    sDist                                    = min(graphData.maxAssumedSdist, 1./graphData.square_sAff(kk1,colorNeighbors) - 1);
    D                                        = exp( -opts_recon.cFactor .* (D{1}(2:end).^2)) .* exp( -graphData.s .* (sDist.^2));
    ssNeighbors                              = find(graphData.square_sAff(kk1,:)>1/graphData.spatialNhoodRadius);
    ssNeighbors                              = setdiff(ssNeighbors, colorNeighbors);
    cDist                                    = pdist2(graphData.colorsForProximal(kk1, :), graphData.colorsForProximal(ssNeighbors, :));
    validssNeighbors                         = find(cDist < graphData.maxColorRadiusForProximal*sqrt(ch/4)); % 20
    if ~isempty(validssNeighbors)
      ssNeighbors                            = ssNeighbors(validssNeighbors);
      cDist                                  = cDist(validssNeighbors);
      colorNeighbors                         = [colorNeighbors ssNeighbors];
      sDist                                  = min(graphData.maxAssumedSdist, 1./graphData.square_sAff(kk1,ssNeighbors) - 1);
      D                                      = [D exp( -opts_recon.cFactor .* (cDist.^2)).*exp( -graphData.s .* (sDist.^2))];
    end
  else
    colorNeighbors                           = find(graphData.square_sAff(kk1,:)>1/graphData.spatialNhoodRadius);
    cDist                                    = pdist2(graphData.colorsForProximal(kk1, :), graphData.colorsForProximal(colorNeighbors, :));
    validcolorNeighbors                      = find(cDist < graphData.maxColorRadiusForProximal*sqrt(ch/4));
    if numel(validcolorNeighbors)<graphData.minEdgeCountForProximal % 2
      [colorNeighbors2, cDist2]              = knnsearch(Mdl, graphData.colorsForDistal(kk1, :), 'K', graphData.minEdgeCountForProximal-numel(validcolorNeighbors));
      colorNeighbors2                        = pure(colorNeighbors2);
      toRemove                               = ismember(colorNeighbors2, colorNeighbors);
      colorNeighbors2(toRemove)              = [];
      cDist2(toRemove)                       = [];
      colorNeighbors                         = [colorNeighbors(validcolorNeighbors) colorNeighbors2];
      cDist                                  = [cDist(validcolorNeighbors) cDist2];
    else
      colorNeighbors                         = colorNeighbors(validcolorNeighbors);
      cDist                                  = cDist(validcolorNeighbors);
    end
    D                                        = [];
    if ~isempty(colorNeighbors)
      sDist                                  = min(graphData.maxAssumedSdist, 1./graphData.square_sAff(kk1,colorNeighbors) - 1);
      D                                      = exp( -opts_recon.cFactor .* (cDist.^2)) .* exp( -graphData.s .* (sDist.^2));
    end
  end
  xx_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = kk1;
  yy_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = colorNeighbors;
  ss_cAff(cIdx:cIdx+numel(colorNeighbors)-1) = D;
  cIdx                                       = cIdx+numel(colorNeighbors);
end

largeSVs                                     = setdiff(1:cc, graphData.smallSVs);
for kk = 1:numel(largeSVs)
  kk1                                        = largeSVs(kk);
  if ccIDsOfSVs(kk1)<=0
    disp('SOMETHING WRONG')
  end
  colorNeighbors                             = setdiff(seedSets{ccIDsOfSVs(kk1)}, kk1, 'stable');
  D                                          = ones(1, numel(colorNeighbors));
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
