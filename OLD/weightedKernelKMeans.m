function index = weightedKernelKMeans(affinityMatrix, index, opts)

% compute the kernel
svCount                  = size(affinityMatrix, 1);
clusterCount             = numel(unique(index));
degreeMatrix             = sparse(1:svCount, 1:svCount, sum(affinityMatrix));
invdegMatrix             = sparse(1:svCount, 1:svCount, 1./(sum(affinityMatrix)+eps));
weightVector             = diag(degreeMatrix);
kernelMatrix             = invdegMatrix * affinityMatrix * invdegMatrix;
index                    = index(:);
oldIndex                 = zeros(size(index));
iter                     = 0;
while any(oldIndex~=index) & iter < opts.maxIter
  DD                     = zeros(svCount, clusterCount);
  for idx = 1:clusterCount
    tempVec              = -2*kernelMatrix(:, index==idx) * weightVector(index==idx);
    tempClusterScalar1   = sum(weightVector(index==idx));
    tempClusterScalar2   = transpose(weightVector(index==idx)) * kernelMatrix(index==idx, index==idx) * weightVector(index==idx);
    DD(:, idx)           = tempVec/tempClusterScalar1 + tempClusterScalar2/(tempClusterScalar1^2);
  end
  DD                     = DD + repmat(diag(kernelMatrix), 1, clusterCount);
  oldIndex               = index;
  [~, index]             = min(DD, [], 2);
  iter                   = iter + 1;
end
