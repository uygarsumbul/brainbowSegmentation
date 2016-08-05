function [tmp, tmpColor] = singleNeuronReconstruction(bbVol, threshold)
% size(bbVol) = [numYpixels,numXpixels,numSlices]; ie ONE CHANNEL
% threshold should be set from 0 to 1

% normalize intensity
bbVol = bbVol - min(bbVol(:));
bbVol = bbVol / max(bbVol(:));

% threshold bbVol
binarybbVol = bbVol > threshold;

%{
getting rid of everything except the largest connected component is a slight gamble,
because it's possible we'll see the single cell split into a few different 
large components.. or one big one and several small ones.. 
also, is background handled here?
%}
CC = bwconncomp(binarybbVol, 26); % get a structure with information about all the connected components in the volume
compsizes = cellfun(@numel,CC.PixelIdxList); % create a list of the sizes of each connected component in the volume
[~,order] = sort(compsizes,'descend'); % find the order of that list, were it to be sorted from biggest to smallest
tmp = false(size(bbVol,1),size(bbVol,2),size(bbVol,3)); % create a blank array the size of the original volume
tmp(CC.PixelIdxList{order(1)}) = true; % mark the indices of all the indices in the largest connected component

iterCount = 0;
while true
    
    neighbors = setdiff(find(imdilate(tmp, true(5,5,5))), find(tmp));
    toRemove = [];
    
    for kk = 2:numel(order)
        if ~isempty(intersect(neighbors, CC.PixelIdxList{order(kk)}))
            tmp(CC.PixelIdxList{order(kk)}) = true;
            toRemove(end+1) = order(kk);
        end
    end
    
    if ~isempty(toRemove)
        CC.PixelIdxList(toRemove) = [];
        compsizes = cellfun(@numel,CC.PixelIdxList);
        [~,order] = sort(compsizes,'descend');
    else
        break;
    end
    iterCount = iterCount+1;
    disp(sprintf('Iteration %d.',iterCount))
end

CCtmp = bwconncomp(tmp, 26);  

myColors   = distinguishable_colors(CCtmp.NumObjects, [0 0 0]); 
tmpColor   = zeros([size(tmp),3]); 
voxelCount = numel(tmp);

for kk=1:CCtmp.NumObjects
    tmpColor(CCtmp.PixelIdxList{kk})                = myColors(kk,1);
    tmpColor(CCtmp.PixelIdxList{kk} + voxelCount)   = myColors(kk,2);
    tmpColor(CCtmp.PixelIdxList{kk} + 2*voxelCount) = myColors(kk,3); 
end


end
