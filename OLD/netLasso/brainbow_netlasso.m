load('netLasso_sv_dcai2_sigma2000_hmin26_0.008_26n_split5e-11_smallCompsRemoved_aff.mat')
load ('svCellsForSuraj_dcai2_sigma2000_hmin26_0.008_26n_split5e-11_smallCompsRemoved.mat');
xlim = stackSize(1)/5; ylim = stackSize(2)/5; zlim = stackSize(3)/5;
size_SVs = {}
chosen_SVs = {};
for i = 1:length(superVoxelCells)
	temp = superVoxelCells(i);
	[a b c] = ind2sub(stackSize, temp{1});
	if any((a<xlim).*(b<ylim).*(c<zlim))
		chosen_SVs{end+1} = i;
		size_SVs{end+1} = length(temp{1});
	end
end
chosen_SVs = cell2mat(chosen_SVs);
size_SVs = cell2mat(size_SVs);size_SVs = size_SVs/sum(size_SVs);
affMat = square_sAff(chosen_SVs, chosen_SVs);
superVoxelMeans = superVoxelMeans(chosen_SVs,:);


[row col val] = find(affMat);
numClosestSVs = 3;
for nodes = 1:size(affMat,1)
	if ~any(nodes==row)
		disp(nodes);
		D = pdist2(superVoxelMeans,superVoxelMeans(nodes,:));
		[B I] = sort(D);
		closestSVindex = I(2:(1+numClosestSVs));
		affMat(nodes,closestSVindex) = 1/numClosestSVs;
	end
end

[row col val] = find(affMat);

M = [row col];
dlmwrite('edges_supervoxel.txt',M,'precision','%i','delimiter',' ');
M = [row col val];
dlmwrite('edges_weights.csv',M,'precision',7,'delimiter',',');
csvwrite('superVoxelMeans.csv',superVoxelMeans);
csvwrite('superVoxelSizes.csv',size_SVs);
