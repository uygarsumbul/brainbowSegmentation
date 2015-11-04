function bigIm = showClusterProjectionsCollage(segmentation, colorStack, options)

if nargin<3; options = []; end;
if ~isfield(options,'displayStyle') || isempty(options.displayStyle); options.displayStyle = 'xy';                 end;
if ~isfield(options,'clusterCount') || isempty(options.clusterCount); options.clusterCount = max(segmentation(:)); end;

% segmentation is a stack with integer values for each cluster
voxelCount = numel(segmentation);
xTileCount = ceil(sqrt(options.clusterCount/40) * 5);
yTileCount = ceil(options.clusterCount/xTileCount);
compx = size(segmentation,1)+size(segmentation,3);
compy = size(segmentation,2)+size(segmentation,3);
switch options.displayStyle
case 'xy'
  bigIm = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, 3);
  for kk=1:options.clusterCount
    [xTile, yTile] = ind2sub([xTileCount yTileCount], kk);
    tmp2 = find(segmentation==kk);
    for mm = 1:3
      tmp1 = zeros(size(segmentation));
      tmp1(tmp2) = colorStack(tmp2+(mm-1)*voxelCount);
      bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, mm) = max(tmp1, [], 3);
    end
  end;
case 'xyz'
  bigIm = ones( (1+compx)*xTileCount-1, (1+compy)*yTileCount-1, 3);
  for kk=1:options.clusterCount;
    tmp1 = colorStack; tmp2 = find(segmentation~=kk); tmp1([tmp2; tmp2+voxelCount; tmp2+2*voxelCount]) = 0;
    tmpxy = squeeze(max(tmp1,[],3));
    tmpzy = permute(squeeze(max(tmp1,[],1)), [2 1 3]);
    tmpxz = squeeze(max(tmp1,[],2));
    thisComponent = ones(compx, compy, 3);
    thisComponent(1:size(segmentation,1), 1:size(segmentation,2), :) = tmpxy;
    thisComponent(size(segmentation,1)+1:end, 1:size(segmentation,2), :) = tmpzy;
    thisComponent(1:size(segmentation,1), size(segmentation,2)+1:end, :) = tmpxz;
    [xTile, yTile] = ind2sub([xTileCount yTileCount], kk);
    bigIm((xTile-1)*(1+compx)+1:xTile*(1+compx)-1, (yTile-1)*(1+compy)+1:yTile*(1+compy)-1, :) = thisComponent;
  end;
case 'xy_meanSubtracted'
  bigIm = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, 3);
  for kk=1:options.clusterCount;
    tmp1 = colorStack; tmp2 = find(segmentation~=kk); tmp1([tmp2; tmp2+voxelCount; tmp2+2*voxelCount]) = 0;
    tmp2 = find(segmentation==kk);
    tmp1(tmp2)              = tmp1(tmp2)              - mean(tmp1(tmp2));
    tmp1(tmp2+voxelCount)   = tmp1(tmp2+voxelCount)   - mean(tmp1(tmp2+voxelCount));
    tmp1(tmp2+2*voxelCount) = tmp1(tmp2+2*voxelCount) - mean(tmp1(tmp2+2*voxelCount));
    tmp = squeeze(max(abs(tmp1),[],3)); % ABSOLUE VALUE
    [xTile, yTile] = ind2sub([xTileCount yTileCount], kk);
    bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, :) = tmp;
  end;
end
%figure; imshow(bigIm,[]);



%bigIm = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, 3);
%for kk=1:clusterCount;
%  voxelCount = numel(segmentation);
%  tmp1 = colorStack; tmp2 = find(segmentation~=kk); tmp1([tmp2; tmp2+voxelCount; tmp2+2*voxelCount]) = 0;
%  tmp = squeeze(max(tmp1,[],3));
%  [xTile, yTile] = ind2sub([xTileCount yTileCount], kk);
%  bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, :) = tmp;
%end;
