function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, vdpgmAssignments, duplications] = demixSupervoxels6(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, opts, origIndex, vdpgmAssignments)
% origIndex must satisfy: setdiff(unique(origIndex), 1:max(origIndex))=[] and setdiff(1:max(origIndex), unique(origIndex))=[]
svMeansNorm                    = svMeans ./ repmat(sqrt(sum(svMeans.^2, 2)), 1, size(svMeans, 2));
allTriplets                    = nchoosek(1:size(svMeans, 2), 3);
allColors                      = zeros(size(svMeans, 1), 3*size(allTriplets, 1));
for kk = 1:size(allTriplets, 1)
  allColors(:, 3*kk-2:3*kk)    = rgb2luv(svMeansNorm(:, allTriplets(kk, :))')';
end
[coeff,score,latent]           = pca(allColors);
svMeansNormLUV                 = score(:, 1:size(svMeans, 2));

binsaff                        = (square_sAff>1/(sqrt(3)+eps));
maxPairResNormSq               = (opts.maxSimilarNeighborNormLUVDist/opts.minImprovementFactor)^2;
duplications                   = zeros(numel(svCells), 3);
voxelCounts                    = voxelCounts(:);
candidates                     = [];
idx                            = 1;
for sv = 1:numel(svCells)
  neighbors                    = find(binsaff(sv, :));
  % IF THERE ARE MORE THAN ONE NEIGHBORS AND THE SUPERVOXEL'S COLOR IS FAR FROM ALL OF THEM
  if numel(neighbors)>1 & voxelCounts(sv)<opts.maxSizeForDemixing & min(pdist2(svMeansNormLUV(sv, :), svMeansNormLUV(neighbors, :))) > opts.maxSimilarNeighborNormLUVDist
    % VALID NEIGHBORS ARE THE NEIGHBORING SUPERVOXELS N SUCH THAT N HAS AT LEAST ONE NEIGHBOR WHOSE COLOR IS SIMILAR TO N'S COLOR 
    validNeighbors             = false(size(neighbors));
    for nn = 1:numel(neighbors)
      thisN                    = neighbors(nn);
      NofN                     = find(binsaff(thisN, :));
      validNeighbors(nn)       = (min(pdist2(svMeansNormLUV(thisN, :), svMeansNormLUV(NofN, :))) < opts.maxSimilarNeighborNormLUVDist);
    end
    neighbors                  = neighbors(validNeighbors);
    % IF MORE THAN ONE VALID (RELIABLE) SOURCES EXIST, TRY DEMIXING
    if numel(neighbors)>1 
      neighborPairs            = nchoosek(neighbors, 2);
      pairResiduals            = zeros(size(neighborPairs, 1), 1);
      for np = 1:size(neighborPairs, 1)
        [~, resnorm, ~]        = lsqnonneg(svMeansNormLUV(neighborPairs(np, :), :)', svMeansNormLUV(sv, :)');
        pairResiduals(np)      = resnorm;
      end
      [mini, pos]              = min(pairResiduals);
      if mini<maxPairResNormSq
        duplications(idx, :)   = [sv neighborPairs(pos, :)];
        idx                    = idx + 1;
      end
    end
  end
end

duplications(idx:end, :)           = []; % disp(size(duplications, 1));
svMeans(duplications(:,1), :)      = [];
for kk = 1:size(duplications, 1)
  svCells{duplications(kk,2)}      = [svCells{duplications(kk,2)}; svCells{duplications(kk,1)}];
  svCells{duplications(kk,3)}      = [svCells{duplications(kk,3)}; svCells{duplications(kk,1)}];
end
voxelCounts(duplications(:,2))     = voxelCounts(duplications(:,2)) + voxelCounts(duplications(:,1));
voxelCounts(duplications(:,3))     = voxelCounts(duplications(:,3)) + voxelCounts(duplications(:,1));
svCells(duplications(:,1))         = [];
voxelCounts(duplications(:,1))     = [];
svColorMins(duplications(:,1), :)  = [];
svColorMaxs(duplications(:,1), :)  = [];
ToRemoveOriginal                   = ismember(origIndex, duplications(:,1));
origIndex(ToRemoveOriginal)        = [];
vdpgmAssignments(ToRemoveOriginal) = [];
mapDown                            = zeros(max(origIndex), 1);
uniqueIndices                      = unique(origIndex);
mapDown(uniqueIndices)             = 1:numel(uniqueIndices);
origIndex                          = mapDown(origIndex);

[row, col, val]                = find(square_sAff);
toUse                          = (ismember(row, duplications(:, 1)) | ismember(col, duplications(:, 1)));
C2                             = [1:size(square_sAff, 1)]';
C3                             = [1:size(square_sAff, 1)]';
C2(duplications(:, 1))         = duplications(:, 2);
C3(duplications(:, 1))         = duplications(:, 3);
row                            = [row; C3(row(toUse)); C2(row(toUse)); C2(row(toUse))];
col                            = [col; C2(col(toUse)); C3(col(toUse)); C2(col(toUse))];
val                            = [val; val(toUse);     val(toUse);     val(toUse)    ];
row(toUse)                     = C3(row(toUse));
col(toUse)                     = C3(col(toUse));
toRemove                       = (row==col);
row(toRemove)                  = [];
col(toRemove)                  = [];
val(toRemove)                  = [];

mapDown                        = zeros(size(square_sAff,1), 1);
uniqueIndices                  = setdiff(1:numel(mapDown), duplications(:,1));
mapDown(uniqueIndices)         = 1:numel(uniqueIndices);
row                            = mapDown(row);
col                            = mapDown(col);
ff                             = size(square_sAff, 1)+1;
id                             = row + col*ff;
[uniqueid, ia, ~]              = unique(id);
newRows                        = row(ia);
newCols                        = col(ia);
newVals                        = zeros(size(newRows));
parfor kk = 1:numel(uniqueid)
  newVals(kk)                  = max(val(id==uniqueid(kk)));
end
square_sAff                    = sparse(newRows, newCols, newVals, numel(uniqueIndices), numel(uniqueIndices));

%tmp                            = sort(duplications(:,1), 'descend');
%for kk = 1:numel(tmp)
%  tt                           = (row>tmp(kk));
%  row(tt)                      = row(tt) - 1;
%  tt                           = (col>tmp(kk));
%  col(tt)                      = col(tt) - 1;
%end
%ff                             = size(square_sAff, 1)+1;
%id                             = row + col*ff;
%[uniqueid, ia, ~]              = unique(id);
%newRows                        = row(ia);
%newCols                        = col(ia);
%newVals                        = zeros(size(newRows));
%parfor kk = 1:numel(uniqueid)
%  newVals(kk)                  = max(val(id==uniqueid(kk)));
%end
%square_sAff                    = sparse(newRows, newCols, newVals, size(square_sAff, 1)-numel(tmp), size(square_sAff, 1)-numel(tmp));
