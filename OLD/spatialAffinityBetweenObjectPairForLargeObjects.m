function aff = spatialAffinityBetweenObjectPairForLargeObjects(xyz1, xyz2, min1, max1, min2, max2, opts)

xyz1( xyz1(:,1)<min2(1)-opts.upperBound, : ) = [];
xyz1( xyz1(:,2)<min2(2)-opts.upperBound, : ) = [];
xyz1( xyz1(:,3)<min2(3)-opts.upperBound, : ) = [];
xyz1( xyz1(:,1)>max2(1)+opts.upperBound, : ) = [];
xyz1( xyz1(:,2)>max2(2)+opts.upperBound, : ) = [];
xyz1( xyz1(:,3)>max2(3)+opts.upperBound, : ) = [];
xyz2( xyz2(:,1)<min1(1)-opts.upperBound, : ) = [];
xyz2( xyz2(:,2)<min1(2)-opts.upperBound, : ) = [];
xyz2( xyz2(:,3)<min1(3)-opts.upperBound, : ) = [];
xyz2( xyz2(:,1)>max1(1)+opts.upperBound, : ) = [];
xyz2( xyz2(:,2)>max1(2)+opts.upperBound, : ) = [];
xyz2( xyz2(:,3)>max1(3)+opts.upperBound, : ) = [];

if isempty(xyz1) | isempty(xyz2)
  aff = 0;
else
  if size(xyz1,1)*size(xyz2,1)<1e6
    aff = 1/min(min(pdist2(xyz1,xyz2)));
  else
    dist = pdist2(xyz1(1,:),xyz2(1,:));
    if size(xyz1,1)>size(xyz2,1)
      for kk = 1:size(xyz2,1)
        dist = min(dist, min(pdist2(xyz1, xyz2(kk, :))));
      end
    else
      for kk = 1:size(xyz1,1)
        dist = min(dist, min(pdist2(xyz1(kk, :), xyz2)));
      end
    end
    aff = 1/dist;
  end
end


