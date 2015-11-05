function im = showIndividualSuperVoxelsXYZ(bbVol, superVoxelCells, superVoxelSubset)

tmpStack = 0*bbVol; 
voxCount = size(bbVol);
voxCount = prod(voxCount(1:3));
compx = size(bbVol,1)+size(bbVol,3)+1;
compy = size(bbVol,2)+size(bbVol,3)+1;
for kk = 1:numel(superVoxelSubset)
  voxelPositions = superVoxelCells{superVoxelSubset(kk)};
  tmpStack([voxelPositions; voxelPositions+voxCount; voxelPositions+2*voxCount]) = bbVol([voxelPositions; voxelPositions+voxCount; voxelPositions+2*voxCount]);
end

tmpxy = squeeze(max(tmpStack,[],3));
tmpzy = permute(squeeze(max(tmpStack,[],1)), [2 1 3]);
tmpxz = squeeze(max(tmpStack,[],2));
thisComponent = ones(compx, compy, 3);
thisComponent(1:size(bbVol,1), 1:size(bbVol,2), :) = tmpxy;
thisComponent(size(bbVol,1)+2:end, 1:size(bbVol,2), :) = tmpzy;
thisComponent(1:size(bbVol,1), size(bbVol,2)+2:end, :) = tmpxz;
im = thisComponent;
