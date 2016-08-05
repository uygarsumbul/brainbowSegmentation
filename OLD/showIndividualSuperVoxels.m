function im = showIndividualSuperVoxels(bbVol, superVoxelCells, superVoxelSubset)

tmpStack = zeros(size(bbVol)); voxCount = size(bbVol); voxCount = prod(voxCount(1:3));
for kk = 1:numel(superVoxelSubset)
  voxelPositions = superVoxelCells{superVoxelSubset(kk)};
  tmpStack([voxelPositions; voxelPositions+voxCount; voxelPositions+2*voxCount]) = bbVol([voxelPositions; voxelPositions+voxCount; voxelPositions+2*voxCount]);
end

im = squeeze(max(tmpStack,[],3));
  
