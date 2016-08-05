function L = assign6Boundaries(L, bbVol)

voxelCount = numel(L);
dims       = size(L);
boundaries = find(L==0);
newBoundaryAssignments = zeros(numel(boundaries), 1);
%% ASSIGN EACH VOXEL ON THE WATERSHED BOUNDARY TO A 6-NEIGHBOR OBJECT (INCLUDING BACKGROUND) BASED ON COLOR
%% WATERSHED DOMAINS MUST BE CALCULATED BY 6-CONNECTIVITY
iterCount = 1;
while ~isempty(boundaries)
  for kk = 1:numel(boundaries)
    thisVox                    = boundaries(kk);
    thisVal                    = [];
    for dd = 1:size(bbVol, 4)
      thisVal                  = [thisVal bbVol(thisVox + (dd-1)*voxelCount)];  
    end
    [xx yy zz]                 = ind2sub(dims, thisVox);
    validNeighbors             = valid26Neighbors(xx, yy, zz, dims);
    validIdx                   = sub2ind(dims, validNeighbors(:,1), validNeighbors(:,2), validNeighbors(:,3));
    validIdx(L(validIdx)==0)   = []; %% BOUNDARY VOXELS ARE INELIGIBLE
    if ~isempty(validIdx)
      validVal                   = [];
      for dd = 1:size(bbVol, 4)
        validVal                 = [validVal bbVol(validIdx + (dd-1)*voxelCount)];
      end
      allNormSq                  = sum((validVal - repmat(thisVal, numel(validIdx), 1)).^2, 2);
      [mini,pos]                 = min(allNormSq);
      newBoundaryAssignments(kk) = L(validIdx(pos));
    end
  end
  L(boundaries) = newBoundaryAssignments;
  boundaries = find(L==0);
  newBoundaryAssignments = zeros(numel(boundaries), 1);
  iterCount = iterCount+1;
end
disp(iterCount)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tmp = valid26Neighbors(xx, yy, zz, dims)
tmp = [[xx-1 yy-1 zz-1]; [xx-1 yy-1 zz]; [xx-1 yy-1 zz+1]; [xx-1 yy zz-1]; [xx-1 yy zz]; [xx-1 yy zz+1]; [xx-1 yy+1 zz-1]; [xx-1 yy+1 zz]; [xx-1 yy+1 zz+1]; ...
       [xx   yy-1 zz-1]; [xx   yy-1 zz]; [xx   yy-1 zz+1]; [xx   yy zz-1];               [xx   yy zz+1]; [xx   yy+1 zz-1]; [xx   yy+1 zz]; [xx   yy+1 zz+1]; ...
       [xx+1 yy-1 zz-1]; [xx+1 yy-1 zz]; [xx+1 yy-1 zz+1]; [xx+1 yy zz-1]; [xx+1 yy zz]; [xx+1 yy zz+1]; [xx+1 yy+1 zz-1]; [xx+1 yy+1 zz]; [xx+1 yy+1 zz+1]];
tmp(:, 1) = min(max(1, tmp(:, 1)), dims(1));
tmp(:, 2) = min(max(1, tmp(:, 2)), dims(2));
tmp(:, 3) = min(max(1, tmp(:, 3)), dims(3));
