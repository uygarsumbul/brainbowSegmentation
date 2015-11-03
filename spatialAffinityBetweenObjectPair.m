function aff = spatialAffinityBetweenObjectPair(xyz1, xyz2, opts)

min1 = min(xyz1, [], 1);
max1 = max(xyz1, [] ,1);
min2 = min(xyz2, [], 1);
max2 = max(xyz2, [], 1);
simpleLowerBound = sqrt(sum(min([min1-min2; min1-max2; max1-min2; max1-max2].^2)));
aff = 0;
if simpleLowerBound<opts.upperBound
  aff = 1/min(min(pdist2(xyz1,xyz2)));
end


  


