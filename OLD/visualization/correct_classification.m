function [svTraced] = correct_classification(bbVol, svTraced, stackSize, svCells, svxyz)
%% svTraced is a list of length number of clusters and each element of list contains a vector of supervoxel IDs
%% svxzy is the xyz coordinates of voxels in each supervoxel. The output of the function: svXYZ.m

%% Returns: updated list of supervoxel clusters svTraced

%% How does it work:
% First it plots the projection collage of current cluster. Click twice on the
% figure at the cluster. First click should be at the cluster (call A) you want to remove the supervoxels from
% and second click should be at the cluster (call B) you want to add the supervoxels
% to. Then, the code will plot the cluster you want to remove the
% supervoxels from. Click on the SVs of cluster A that you want to assign
% to the cluster B(Multiple clicks are allowed). Then, the code will
% do the reassignment of the SVs. It'll then ask if you want to continue.
% Press 0 if you don't and 1 if you do. It'll repeat that until you press
% 0. The updated cluster of supervoxels are returned. 



%% Code
bbVol=bbVol/max(bbVol(:));

mySeg = zeros(stackSize);
for kk1=1:length(svTraced) %
	if(numel(svTraced{kk1})>0)
		for kk2=1:numel(svTraced{kk1})
			mySeg(svCells{svTraced{kk1}(kk2)})=kk1;
		end
	end
end

clusterCount = length(svTraced);
xTileCount = ceil(sqrt(clusterCount/40) * 5);

more = 1;
while more
    bigIm = showClusterProjectionsCollage(mySeg, bbVol);
    figure(1);imshow(bigIm,[]);
	[y,x] = ginput(2);
    close; 
	x_tile_from = ceil(x(1)/stackSize(1)); % x increases as we go down the image
	y_tile_from = ceil(y(1)/stackSize(2));
	x_tile_to = ceil(x(2)/stackSize(1));
	y_tile_to = ceil(y(2)/stackSize(2));
	from_cluster_num = x_tile_from+(y_tile_from-1)*xTileCount;
	to_cluster_num = x_tile_to+(y_tile_to-1)*xTileCount;
	single_clus_plot = figure(2);
	im = showIndividualSuperVoxelsXYZ(bbVol,svCells,svTraced{from_cluster_num}); figure(2);imshow(im,[]);
	[y,x] = ginput();
    close;
	point = round([x,y]);
	goodSV = {};
    for iter=1:length(x)
        % should we do distance or just matching?
        distances = zeros(length(svTraced{from_cluster_num}),1);
        for id = 1:length(distances)
            distances(id)=sqrt(sum((point(iter,:)-svxyz(svTraced{from_cluster_num}(id),1:2)).^2));
        end
        found = 0;
        [sortD indSort] = sort(distances);
        temp = 0;
        while ((found==0) && (temp<length(svTraced{from_cluster_num})))
            closeSV = svTraced{from_cluster_num}(indSort(temp+1));
            closeSVvoxels = svCells{closeSV};
            for iii = 1:length(closeSVvoxels)
                [x y z]=ind2sub(stackSize,closeSVvoxels(iii));
                if (sum(point(iter,:)==[x,y])==2)
                    found=found+1;
                    disp('found'); disp(closeSV)
                    goodSV{length(goodSV)+1}=closeSV;
                    break;
                end
            end
            temp = temp+1;
        end
    end
    for iter = 1:length(goodSV)
    	svTraced{from_cluster_num} = svTraced{from_cluster_num}(svTraced{from_cluster_num}~=goodSV{iter});
    	svTraced{to_cluster_num}(numel(svTraced{to_cluster_num})+1) = goodSV{iter};
    end
    
    mySeg = zeros(stackSize);
    for kk1=1:length(svTraced)
        if(numel(svTraced{kk1})>0)
            for kk2=1:numel(svTraced{kk1})
                mySeg(svCells{svTraced{kk1}(kk2)})=kk1;
            end
        end
    end
    more = input('Want more?');
end

end