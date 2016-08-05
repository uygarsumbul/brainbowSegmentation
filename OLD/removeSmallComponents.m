function mask = removeSmallComponents(mask, opts)

CC = bwconncomp(mask, 26);
componentSizes = zeros(CC.NumObjects, 1);
parfor kk = 1:CC.NumObjects
  componentSizes(kk) = numel(CC.PixelIdxList{kk});
end
smallComponents = find(componentSizes<opts.minVoxelCount);
toRemove        = false(numel(smallComponents), 1);
parfor kk = 1:numel(smallComponents)
  [xx, yy, zz]       = ind2sub(size(mask), CC.PixelIdxList{smallComponents(kk)});
  xyz                = [xx-mean(xx), yy-mean(yy), zz-mean(zz)];
  zz                 = zz * opts.zAnisotropy; % zAnisotropy = zVoxelLength/xyVoxelLength
  %% FIND THE INERTIA TENSOR
  weights            = ones(size(xyz,1),1);
  inertiaTensor      = zeros(3);
  inertiaTensor(1,1) = sum(weights .* (xyz(:,2).^2 + xyz(:,3).^2)); inertiaTensor(2,2) = sum(weights .* (xyz(:,1).^2 + xyz(:,3).^2));
  inertiaTensor(3,3) = sum(weights .* (xyz(:,1).^2 + xyz(:,2).^2)); inertiaTensor(1,2) = -sum(weights .* xyz(:,1) .* xyz(:,2));
  inertiaTensor(1,3) = -sum(weights .* xyz(:,1) .* xyz(:,3)); inertiaTensor(2,3) = -sum(weights .* xyz(:,2) .* xyz(:,3));
  inertiaTensor(2,1) = inertiaTensor(1,2); inertiaTensor(3,1) = inertiaTensor(1,3); inertiaTensor(3,2) = inertiaTensor(2,3);
  %% FIND THE MOMENT OF INERTIA FOR EACH PRINCIPAL AXIS
  [principalAxes, evMatrix] = eig(inertiaTensor);
  %% IF THE RATIO OF THE LARGEST TO SMALLEST EIGENVALUES IS SMALL, THEN THE SUPERVOXEL IS NOT ELONGATED ENOUGH
  %% IT MAY NOT BE A DENDRITE (MAY BE AN ORGANEL) -- mCherry information can be useful here.
  if evMatrix(3,3)/evMatrix(1,1)<opts.moiRatioThreshold
    toRemove(kk) = true;
  end
end
toRemove = smallComponents(toRemove);
%% REMOVE
for kk = 1:numel(toRemove)
  mask(CC.PixelIdxList{toRemove(kk)}) = 0;
end
