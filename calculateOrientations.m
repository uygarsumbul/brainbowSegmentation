function [orientations, conditionNumbers, meanPos] = calculateOrientations(svCells, stackSize, opts)

if ~isfield(opts,'condStyle') || isempty(opts.condStyle); opts.condStyle = 1; end;

orientations                = zeros(3, numel(svCells));
conditionNumbers            = zeros(1, numel(svCells));
meanPos                     = zeros(3, numel(svCells));
parfor kk = 1:numel(svCells)
  [xx, yy, zz]              = ind2sub(stackSize, svCells{kk});
  zz                        = zz(:) * opts.zAnisotropy; % zAnisotropy = zVoxelLength/xyVoxelLength
  meanPos(:, kk)            = mean([xx, yy, zz], 1)';
  xyz                       = [xx(:)-mean(xx), yy(:)-mean(yy), zz(:)-mean(zz)];
  %% FIND THE INERTIA TENSOR
  weights                   = ones(size(xyz,1),1);
  inertiaTensor             = zeros(3);
  inertiaTensor(1,1)        = sum(weights .* (xyz(:,2).^2 + xyz(:,3).^2));
  inertiaTensor(2,2)        = sum(weights .* (xyz(:,1).^2 + xyz(:,3).^2));
  inertiaTensor(3,3)        = sum(weights .* (xyz(:,1).^2 + xyz(:,2).^2));
  inertiaTensor(1,2)        = -sum(weights .* xyz(:,1) .* xyz(:,2));
  inertiaTensor(1,3)        = -sum(weights .* xyz(:,1) .* xyz(:,3));
  inertiaTensor(2,3)        = -sum(weights .* xyz(:,2) .* xyz(:,3));
  inertiaTensor(2,1)        = inertiaTensor(1,2);
  inertiaTensor(3,1)        = inertiaTensor(1,3);
  inertiaTensor(3,2)        = inertiaTensor(2,3);
  %% FIND THE MOMENT OF INERTIA FOR EACH PRINCIPAL AXIS
  [principalAxes, evMatrix] = eig(inertiaTensor);
  evVal                     = diag(evMatrix);
  evVal(evVal<1e-5)         = 1e-5;
  orientations(:, kk)       = principalAxes(:, 1);
  if opts.condStyle == 1
    conditionNumbers(kk)    = evVal(3)/evVal(1);
  else
    conditionNumbers(kk)    = evVal(3)/evVal(2);
  end
end
