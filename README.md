# Welcome
## Overview
* A presentation on the basics of the approach is available here: http://www.stat.columbia.edu/~uygar/dapg2015verbose.pdf

* denoisingScript: denoise brainbow stack

* superVoxelizationScript: switch to an over segmented 'supervoxel' representation

* colorKmeans: k-means clustering based on color (LUV space)

* calculateNormalizedCuts_mem: semi-supervised normalized cuts segmentation of the supervoxels

* manipulateClusters: cluster-level operations (human interaction) to provide cluster seeds for semi-supervised graph cuts i.e., merge, split, clean

* retainLargeClumps: retain supervoxel subsets only if they form large enough connected components (human interaction) to provide cluster seeds for semi-supervised graph cuts.

* visualization: directory including visualization and file-write scripts

    - segProjectionWrite_script: script to visualize maximum intensity z-projections of all clusters in one big 2d image.

    - segProjectionWrite_scriptXYZ: script to visualize maximum intensity x-, y-, and z-projections of all clusters in one big 2d image.

    - write3dTifForIndividualComponents_script: script to generate a 3d tif file where only voxels belonging to a particular component are shown

## Help
* Denoising is performed by the BM4D collaborative filtering package, which is available here:
http://www.cs.tut.fi/~foi/GCF-BM3D/index.html#ref_software

## To-do list
* SLIC is another way of supervoxelization whose supervoxels are more uniform in size compared to the current watershed based approach.
The code for 2d super pixels is available here: http://www.mathworks.com/matlabcentral/answers/uploaded_files/18687/slic.m
* Network LASSO can be used to remove local color variations. The code is available here:
http://snap.stanford.edu/snapvx/
* Other ways of graph partitioning such as affinity propagation can be used instead of normalized cuts.
* Better visualization options/scripts within MATLAB

## How/What to run
* Step 0: denoise the raw image stack using denoisingScript.m. sigma=2000 is an overestimate of noise standard deviation that seems to perform well.
bbVol is a 4-d array of size xLength x yLength x zLength x channelCount.
Simple parallelization over individual channels can be achieved by changing 'for' into 'parfor' when memory is not an issue. (~ a few hours for a 1000x1000x200 stack)
* Step 1: identify the foreground and generate a supervoxel representation using superVoxelizationScript.m.
Various options (sensitivity in identifying foreground, spatial distance calculation options,...) are set at the beginning of the script.
If memory is not an issue, some steps can be parallelized by simply changing 'for's into 'parfor's. (~6 hours for a 1000x1000x200 stack)
* Step 2: perform a basic k-means clustering based on color information only. Typically, results improve significantly if colors are represented in the LUV color space.
Here, the variables svMeans and svCells are computed in Step 1. clusterCount indicates (an estimate of) the number of neurons in the stack. The output index is an array indicating
the cluster each supervoxel belongs to. (~a few minutes for a 1000x1000x200 stack)
The kmeans algorithm can be parallelized as described here: http://theory.stanford.edu/~sergei/papers/vldb12-kmpar.pdf

cc                                           = size(svMeans, 1);

voxCounts                                    = zeros(cc, 1);

for kk = 1:cc

  voxCounts(kk)                              = numel(svCells{kk});

end

opts_fkmeans.weight                          = sqrt(voxCounts);

opts_fkmeans.careful                         = true;

opts_fkmeans.maxIter                         = 200;

if size(svMeans,2)==3

  svMeansLUV                                 = rgb2luv(svMeans')';

else

  svMeansLUV(:,1:3)                          = rgb2luv(svMeans(:,1:3)')';

  svMeansLUV(:,4:end)                        = svMeans(:,4:end)*50;

end

index                                        = colorKmeans(clusterCount, svMeansLUV, opts_fkmeans);

* Step 3: Cluster-level manipulation of the initial results. The function manipulateClusters.m can be used multiple times. retainLargeClumps.m removes small, disconnected
supervoxels in individual clusters because they are more likely to represent mistakes. The remaining 'cleaned-up' clusters can be used as seed clusters in further clustering.

** Example: remove cluster 5, split clusters 1 and 9 into 2, merge clusters 15 and 16, and merge clusters 13 and 18.

---manipulationSets.removalSet = [5]; manipulationSets.splitSet = [1 9]; manipulationSets.splitCount = 2; manipulationSets.mergeSets = {[15 16], [13 18]}; 
---index = snrAwareUserInteraction(index, manipulationSets, graphData.colorData, opts_fkmeans);
