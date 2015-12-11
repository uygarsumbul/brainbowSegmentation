function index = nCutsForSCMerging(binSaff, colorData)

cc                                       = size(colorData, 1);
if cc == 2
  index                                  = [1 2];
  return;
end

opts_fkmeans.careful                     = true;
opts_fkmeans.maxIter                     = 1000;
opts_irbleigs.K                          = 2;
opts_recon.cDecay                        = 1e-2;
index                                    = zeros(cc, 1);
xx_cAff                                  = zeros(round(0.1*cc*(cc-1)/2), 1);
yy_cAff                                  = zeros(round(0.1*cc*(cc-1)/2), 1);
ss_cAff                                  = zeros(round(0.1*cc*(cc-1)/2), 1);
cIdx                                     = 1;
%% LOOP OVER INDIVIDUAL SUPERVOXELS
for kk1 = 1:cc
  sNeighbors                             = find(binSaff(kk1, 1:kk1-1));
  xx_cAff(cIdx:cIdx+numel(sNeighbors)-1) = kk1;
  yy_cAff(cIdx:cIdx+numel(sNeighbors)-1) = sNeighbors;
  ss_cAff(cIdx:cIdx+numel(sNeighbors)-1) = exp(-opts_recon.cDecay .* (pdist2(colorData(kk1,:), colorData(sNeighbors, :)).^2));
  cIdx                                   = cIdx+numel(sNeighbors);
end
unused                                   = find(~xx_cAff, 1, 'first');
if ~isempty(unused)
  ss_cAff(unused:end)                    = [];
  xx_cAff(unused:end)                    = [];
  yy_cAff(unused:end)                    = [];
end
affinityMatrix                           = sparse(xx_cAff, yy_cAff, ss_cAff, cc, cc);
affinityMatrix                           = max(affinityMatrix, transpose(affinityMatrix));
%% CALCULATE TOP 2 EIGENVECTORS OF THE LAPLACIAN
DD                                       = sparse(1:cc, 1:cc, 1./sqrt(sum(affinityMatrix)+eps));
if cc<1000
  try
    [VV, DD]                             = eig(full(DD*affinityMatrix*DD));
    [DDdiag, ii]                         = sort(diag(DD), 'descend');
    top2Eigenvectors                     = VV(:, ii(1:2)); % [top2Eigenvectors, ~] = eigs(DD*affinityMatrix*DD, 2);
  catch
    disp(affinityMatrix)
    disp(size(affinityMatrix))
    index                                = ones(1, cc);
    return;
  end
else
  [top2Eigenvectors, ss, PRGINF]         = irbleigs(DD*affinityMatrix*DD, opts_irbleigs);
end
%% NORMALIZE INDIVIDUAL EIGENVECTORS
for ii=1:size(top2Eigenvectors,1)
  top2Eigenvectors(ii,:)                 = top2Eigenvectors(ii,:) / norm(top2Eigenvectors(ii,:));
end
distortion                               = 1e60;
for kk = 1:200
  [myIndex, myCentroids, myDisto]        = simple_fkmeans(top2Eigenvectors, opts_irbleigs.K, opts_fkmeans);
  if sum(myDisto)<distortion
    distortion                           = sum(myDisto);
    index                                = myIndex;
  end
end

