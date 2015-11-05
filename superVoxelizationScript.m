superVoxelOpts.HMINTH26                                  = 0.008;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound = 31;
superVoxelOpts.colorDistanceUpperBound                   = 0.03;
superVoxelOpts.splitHardSVopts.detThreshold              = 5e-11;
superVoxelOpts.splitHardSVopts.connectivity              = 26;
superVoxelOpts.splitHardSVopts.subdivisionSizeThreshold  = 50;
superVoxelOpts.removeSmallComponents.minVoxelCount       = 20; %100;
superVoxelOpts.removeSmallComponents.moiRatioThreshold   = 3;
superVoxelOpts.removeSmallComponents.zAnisotropy         = 3;
% superVoxelOpts.processLoneSuperVoxels.colorSim           = 5;   % obsolete
% superVoxelOpts.processLoneSuperVoxels.larsThreshold      = 0.2; % obsolete
superVoxelOpts.dataset                                   = '~/bb/data/dcai_brainbow2_bm4d_sigma2000.mat';
superVoxelOpts.filePreamble                              = 'sv_dcai2_sigma2000_hmin26_0.008_26n_split5e-11_smallCompsRemoved';
load(superVoxelOpts.dataset); bbVol(bbVol<0)=0;
for kk = 1:size(bbVol, 4)
  rawStack = bbVol(:,:,:,kk); rawStack = rawStack - min(rawStack(:)); rawStack = rawStack / max(rawStack(:)); bbVol(:,:,:,kk) = rawStack;
end
clear rawStack;
gradAmplitude = zeros(size(bbVol,1), size(bbVol,2), size(bbVol,3));
tmpbb = zeros(size(bbVol,1)+1,size(bbVol,2),size(bbVol,3),size(bbVol,4)); tmpbb(1:end-1,:,:,:) = bbVol; tmpbb(2:end,:,:,:) = bbVol;
gradAmplitude = max(gradAmplitude, squeeze(max(abs(diff(tmpbb,1,1)),[],4)) );       
tmpbb = zeros(size(bbVol,1),size(bbVol,2)+1,size(bbVol,3),size(bbVol,4)); tmpbb(:,1:end-1,:,:) = bbVol; tmpbb(:,2:end,:,:) = bbVol;
gradAmplitude = max(gradAmplitude, squeeze(max(abs(diff(tmpbb,1,2)),[],4)) );
tmpbb = zeros(size(bbVol,1),size(bbVol,2),size(bbVol,3)+1,size(bbVol,4)); tmpbb(:,:,1:end-1,:) = bbVol; tmpbb(:,:,2:end,:) = bbVol;
gradAmplitude = max(gradAmplitude, squeeze(max(abs(diff(tmpbb,1,3)),[],4)) );
clear tmpbb;
disp('DETECTING AND REMOVING THE BACKGROUND');
%% 26-WATERSHED SEGMENTATION, 0 DENOTES WATERSHED BOUNDARIES, 1 DENOTES BACKGROUND OBJECT
tic; L = watershed(imhmin(gradAmplitude, superVoxelOpts.HMINTH26), 26); toc;
clear gradAmplitude;
stackSize = size(L); voxelCount = prod(stackSize);
backgroundVox = find(L==1);
disp('ASSIGNING WATERSHED BOUNDARIES TO OBJECTS');
tic; L = assign6Boundaries(L, bbVol); toc;
disp('REMOVING SMALL COMPONENTS');
tic; mask = removeSmallComponents(L~=1, superVoxelOpts.removeSmallComponents);
L(~mask) = 1; % 1 IS ASSIGNED TO BACKGROUND IN THIS DATASET - SHOULD REVISE FOR GENERAL CASE
for kk = max(L(:)):-1:2; if ~any(L(:)==kk); L(L>kk) = L(L>kk) - 1; end; end; toc;
disp('FORMING SUPERVOXEL CELLS');
superVoxelCells = cell(1,max(L(:))-1); % 1 IS ASSIGNED TO BACKGROUND IN THIS DATASET - SHOULD REVISE FOR GENERAL CASE
tic; parfor kk = 1:max(L(:))-1; superVoxelCells{kk} = find(L==kk+1); end;  toc;
disp(numel(superVoxelCells))
%% SPLIT HETEROGENEOUS SUPERVOXELS USING COLOR AND CONNECTIVITY - 1 WILL NOT SIGNIFY BACKGROUND ANY MORE
disp('SPLITTING HETEROGENEOUS SUPERVOXELS');
tic; [L, superVoxelCells] = splitHardSuperVoxels(superVoxelOpts.splitHardSVopts, superVoxelCells, bbVol); toc;
disp(numel(superVoxelCells))
disp('CALCULATING SUPERVOXEL MEANS AND BOUNDARIES');
cc                         = numel(superVoxelCells);
superVoxelMeans            = zeros(cc, size(bbVol, 4));
boundaryVoxels             = cell(1, cc);
boundaryVoxelsSub          = cell(1, cc);
tic;
parfor kk = 1:cc
  [xx,yy,zz]               = ind2sub(stackSize, superVoxelCells{kk});
  xSub                     = min(xx)-1;
  ySub                     = min(yy)-1;
  zSub                     = min(zz)-1;
  xx                       = xx-xSub;
  yy                       = yy-ySub;
  zz                       = zz-zSub;
  maxxx                    = max(xx);
  maxyy                    = max(yy);
  maxzz                    = max(zz);
  tmp                      = false(maxxx, maxyy, maxzz);
  reducedIndices           = sub2ind([maxxx, maxyy, maxzz], xx, yy, zz);
  tmp(reducedIndices)      = true;
  localBoundaryVoxels      = find(tmp & ~imerode(tmp, ones(3,3,3)));
  if ~isempty(localBoundaryVoxels)
    [xx,yy,zz]             = ind2sub(size(tmp), localBoundaryVoxels);
  end
  xx                       = xx+xSub;
  yy                       = yy+ySub;
  zz                       = zz+zSub;
  boundaryVoxelsSub{kk}    = [xx,yy,zz];
  boundaryVoxels{kk}       = sub2ind(stackSize, xx, yy, zz);
end
for kk = 1:cc
  thisSVcolors             = [];
  for dd = 1:size(bbVol, 4)
    thisSVcolors           = [thisSVcolors bbVol(superVoxelCells{kk} + (dd-1)*voxelCount)];
  end
  superVoxelMeans(kk,:)    = mean(thisSVcolors, 1);
end
toc;
disp('SAVING SUPERVOXEL BOUNDARIES');
cd ~/bb/data;
fileName = [superVoxelOpts.filePreamble '.mat'];
save(fileName, 'superVoxelOpts', 'superVoxelCells', 'superVoxelMeans', 'boundaryVoxels', 'boundaryVoxelsSub');
cd ~/bb;
disp('CALCULATING SUPERVOXEL AFFINITIES IN SPACE');
tic; sAff = calculate_sAff(cc, boundaryVoxelsSub, superVoxelOpts.spatialDistanceCalculationOpts); toc;
% disp('CALCULATING SUPERVOXEL AFFINITIES IN COLOR');                                           % obsolete
% tic; cAff = calculate_cAff(cc, superVoxelMeans, superVoxelOpts.colorDistanceUpperBound); toc; % obsolete
disp('SAVING');
cd ~/bb/data;
fileName = [superVoxelOpts.filePreamble '_aff.mat'];
save(fileName, 'superVoxelOpts', 'superVoxelCells', 'superVoxelMeans', 'stackSize', 'sAff', 'cAff', 'boundaryVoxels', '-v7.3');
cd ~/bb;

%disp('SPLITTING SUPERVOXELS AT BRANCH INTERSECTIONS')
%tic; [svCells,svMeans,svMeansInt,boundaryVoxels,boundaryVoxelsSub,duplicateSVs] = ...
%    processLoneSuperVoxels(superVoxelCells,superVoxelMeans,superVoxelMeansInt,stackSize,sAff,cAff,boundaryVoxels,boundaryVoxelsSub,superVoxelOpts.processLoneSuperVoxels); toc;
%cc = numel(svCells);
%disp('CALCULATING SUPERVOXEL AFFINITIES IN SPACE AFTER INTERSECTION SPLITTING -- INEFFICIENT');
%tic; sAff = calculate_sAff(cc, boundaryVoxelsSub, superVoxelOpts.spatialDistanceCalculationOpts); toc;
%disp('ASSIGNING ZERO SPATIAL AFFINITY BETWEEN DUPLICATE SUPERVOXELS -- CARE NEEDED IN CLUSTER ASSIGNMENTS');
%for kk = 1:size(duplicateSVs, 1)
%  tmpmin    = min(duplicateSVs(kk, 1:2));
%  tmpmax    = max(duplicateSVs(kk, 1:2));
%  pos       = cc*tmpmin - cc - tmpmin*(tmpmin+1)/2 + tmpmax;
%  sAff(pos) = 0;
%end
%disp('CALCULATING SUPERVOXEL AFFINITIES IN COLOR AFTER INTERSECTION SPLITTING');
%tic; cAff = calculate_cAff(cc, svMeansInt, superVoxelOpts); toc;
%disp('SAVING AFTER INTERSECTIONS ARE SPLIT');
%cd ~/bb/data;
%fileName = ['splitIntersection_' superVoxelOpts.filePreamble '_aff.mat'];
%save(fileName, 'superVoxelOpts', 'svCells', 'svMeans', 'svMeansInt', 'stackSize', 'sAff', 'cAff', 'boundaryVoxels', 'duplicateSVs', '-v7.3');
%cd ~/bb;
