
componentID = 1;
segmentation = false(stackSize);
thisCluster = find(index==componentID);
for kk=1:numel(thisCluster)
  segmentation(svCells{thisCluster(kk)}) = true;
end
theseVoxels = find(segmentation);
tmp = zeros(size(bbVol));
for kk = 1:size(bbVol,4)
  tmp(theseVoxels+(kk-1)*voxelCount) = bbVol(theseVoxels+(kk-1)*voxelCount);
end
imwrite(squeeze(tmp(:, :, 1, :)), ['tif3d_component' num2str(componentID) '.tif']);
for kk = 2:size(bbVol, 3)
  imwrite(squeeze(tmp(:, :, kk, :)), ['tif3d_component' num2str(componentID) '.tif'], 'WriteMode', 'append');
end
