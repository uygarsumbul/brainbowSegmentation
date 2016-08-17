function supervoxelize(superVoxelOpts)
tic;disp('BASIC WATERSHED');
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0;
for kk = 1:size(bbVol, 4)
  rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack;
end
clear rawStack;
channelCount                                                             = size(bbVol, 4);
superVoxelOpts.brightnessThreshold                                       = superVoxelOpts.brightnessThreshold  * sqrt(channelCount/3);
gradAmplitude                                                            = zeros(size(bbVol,1), size(bbVol,2), size(bbVol,3));
gradAmplitude                                                            = max(gradAmplitude, squeeze(max(abs(cat(1, zeros(1,size(bbVol,2),size(bbVol,3),size(bbVol,4)), diff(bbVol,1,1))),[],4)) );
gradAmplitude                                                            = max(gradAmplitude, squeeze(max(abs(cat(2, zeros(size(bbVol,1),1,size(bbVol,3),size(bbVol,4)), diff(bbVol,1,2))),[],4)) );
gradAmplitude                                                            = max(gradAmplitude, squeeze(max(abs(cat(3, zeros(size(bbVol,1),size(bbVol,2),1,size(bbVol,4)), diff(bbVol,1,3))),[],4)) );
L                                                                        = watershed(imhmin(gradAmplitude, superVoxelOpts.HMINTH26), 26); clear gradAmplitude;
[maxi26, pos26]                                                          = max(hist(L(L~=0),1:double(max(L(:))))); bgvoxels26 = find(L==pos26); L(L==1) = pos26; L(bgvoxels26) = 1; % 1 IS THE LARGEST COMPONENT (BACKGROUND)
toc;
tic;disp('ASSIGN BOUNDARIES');
L                                                                        = assign6Boundaries(L, bbVol);
stackSize                                                                = size(L);
voxelCount                                                               = prod(stackSize);
toc;
tic;disp('SUPERVOXEL CELLS - 0 IS THE BACKGROUND, EVERYTHING ELSE IS A REGULAR SUPERVOXEL');
Lmap                                                                     = 0:max(L(:))-1;
L(L>0)                                                                   = Lmap(L(L>0));
[pos, ~, id]                                                             = find(L(:));
[sortedID, idx]                                                          = sort(id);
sortedPos                                                                = pos(idx);
svCount                                                                  = max(L(:));
superVoxelCells                                                          = cell(1, svCount);
for sv = 1:svCount-1
  superVoxelCells{sv}                                                    = sortedPos(find(sortedID==sv, 1):find(sortedID==sv+1, 1)-1);
end
superVoxelCells{end}                                                     = sortedPos(find(sortedID==svCount, 1):end);
toc;
tic;disp('SPLIT HARD SUPERVOXELS');
superVoxelCells                                                          = splitInconsistentSuperVoxels(superVoxelOpts.splitInconsistentSVopts, superVoxelCells, bbVol);
toc; disp(numel(superVoxelCells)); % L NEEDS UPDATING IF YOU WANT TO USE IT FURTHER BEYOND FOREGROUND/BACKGROUND DETECTION
tic;disp('THRESHOLD AND WARP THE THRESHOLDED IMAGE INTO WATERSHED IMAGE')
voxelNorms                                                               = sqrt(sum(bbVol.^2,4));
proposal                                                                 = (voxelNorms>superVoxelOpts.brightnessThreshold);
proposal                                                                 = proposal | (L>0);
voxelNorms(~proposal)                                                    = 0;
voxelNorms                                                               = sparse(voxelNorms(:));
proposal                                                                 = topologyPreservingVolumeShrinkerWithPreference6(proposal, L>0, 30, voxelNorms);
toc;
tic;disp('GENERATE SUPERVOXELS FROM THE DIFFERENCE VOXELS AFTER WARPING AND ADD TO THE EXISTING LIST');
CC                                                                       = bwconncomp(proposal & L==0, 26);
additionalSupervoxels                                                    = cell(1, CC.NumObjects);
for kk = 1:CC.NumObjects
  additionalSupervoxels{kk}                                              = CC.PixelIdxList{kk};
end
additionalSupervoxels                                                    = splitInconsistentSuperVoxels(superVoxelOpts.splitInconsistentSVopts, additionalSupervoxels, bbVol);
superVoxelCells                                                          = [superVoxelCells additionalSupervoxels];
svCount                                                                  = numel(superVoxelCells);
toc; disp(svCount);
tic;disp('CALCULATE BOUNDARY VOXELS');
boundaryVoxels                                                           = cell(1, svCount);
boundaryVoxelsSub                                                        = cell(1, svCount);
parfor kk = 1:svCount
  [xx,yy,zz]                                                             = ind2sub(stackSize, superVoxelCells{kk});
  xSub                                                                   = min(xx)-2;
  ySub                                                                   = min(yy)-2;
  zSub                                                                   = min(zz)-2;
  xx                                                                     = xx-xSub;
  yy                                                                     = yy-ySub;
  zz                                                                     = zz-zSub;
  maxxx                                                                  = max(xx);
  maxyy                                                                  = max(yy);
  maxzz                                                                  = max(zz);
  tmp                                                                    = false(maxxx+1, maxyy+1, maxzz+1);
  reducedIndices                                                         = sub2ind([maxxx+1, maxyy+1, maxzz+1], xx, yy, zz);
  tmp(reducedIndices)                                                    = true;
  localBoundaryVoxels                                                    = find(tmp & ~imerode(tmp, ones(3,3,3)));
  if ~isempty(localBoundaryVoxels)
    [xx,yy,zz]                                                           = ind2sub(size(tmp), localBoundaryVoxels);
  end
  xx                                                                     = xx+xSub;
  yy                                                                     = yy+ySub;
  zz                                                                     = zz+zSub;
  boundaryVoxelsSub{kk}                                                  = [xx,yy,zz];
  boundaryVoxels{kk}                                                     = sub2ind(stackSize, xx, yy, zz);
end
tic;disp('CALCULATE SPATIAL AFFINITIES');
sAff                                                                     = calculate_sAff(svCount, boundaryVoxelsSub, superVoxelOpts.spatialDistanceCalculationOpts);
toc;
tic;disp('CALCULATE SUPERVOXEL MEANS')
superVoxelMeans                                                          = zeros(svCount, channelCount);
for kk = 1:svCount
  thisSVcolors                                                           = zeros(numel(superVoxelCells{kk}), channelCount);
  for dd = 1:channelCount
    thisSVcolors(:, dd)                                                  = bbVol(superVoxelCells{kk} + (dd-1)*voxelCount);
  end
  superVoxelMeans(kk,:)                                                  = mean(thisSVcolors, 1); % superVoxelMeans DOES NOT USE AUGMENTED VOXELS!
end
toc;
tic;disp('CALCULATE SUPERVOXEL COLOR EXTREMES')
svColorMins                                                              = zeros(numel(superVoxelCells), size(bbVol, 4));
svColorMaxs                                                              = zeros(numel(superVoxelCells), size(bbVol, 4));
normalizer                                                               = sqrt(sum(bbVol.^2, 4));
for dd = 1:size(bbVol, 4)
  bbVol(:,:,:,dd) = bbVol(:,:,:,dd) ./ normalizer;
end; clear normalizer;
for kk = 1:numel(superVoxelCells)
  tmp                                                                    = zeros(numel(superVoxelCells{kk}), size(bbVol, 4));
  for dd = 1:size(bbVol, 4)
    tmp(:, dd)                                                           = bbVol(superVoxelCells{kk} + (dd-1)*voxelCount);
  end
  svColorMins(kk, :)                                                     = min(tmp, [], 1);
  svColorMaxs(kk, :)                                                     = max(tmp, [], 1);
end
toc;
tic;disp('SAVE VARIABLES');
fileName                                                                 = [superVoxelOpts.filePreamble '_aff.mat'];
save(fileName, 'superVoxelOpts', 'superVoxelCells', 'superVoxelMeans', 'stackSize', 'sAff', 'boundaryVoxels', 'svColorMins', 'svColorMaxs', '-v7.3');
cd ~/bb;
toc;
