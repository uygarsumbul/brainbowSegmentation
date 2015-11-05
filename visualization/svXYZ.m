% Calculate svXYZ for the list of SVs
function svxyz = svXYZ(svCells, stackSize)

svxyz = zeros(length(svCells),3);
for iij = 1:length(svCells)
    currSVvoxels = svCells{(iij)};
    sumxyz = zeros(1,3);
    for iii = 1:length(currSVvoxels)
        [x y z]=ind2sub(stackSize,currSVvoxels(iii));
        sumxyz = sumxyz+[x y z];
    end
    svxyz((iij),:) = sumxyz/length(currSVvoxels);
end

end

