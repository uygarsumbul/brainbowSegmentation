function groundTruth = svPartitionFromTraces(traceFileNames, svCells, stackSize)

svSeg                                                   = zeros(stackSize);
for kk = 1:numel(svCells)
  svSeg(svCells{kk})                                    = kk;
end
groundTruth                                             = zeros(numel(svCells), 1);
for kk = 1:numel(traceFileNames)
  [nodes, edges, radii, nodeTypes, abort]               = readArborTrace(traceFileNames{kk}, -1:6);
  nodes(:, 1:2)                                         = nodes(:, 2:-1:1) + 1;
  ind                                                   = sub2ind(stackSize, nodes(:,1), nodes(:,2), nodes(:,3));
  groundTruth(setdiff(unique(svSeg(ind)), 0, 'stable')) = kk;
end

% im = showIndividualSuperVoxels(bbVol, svCells, find(groundTruth)); figure; imshow(im(:,:,1:3),[]);

%  unassigned = find(ind==0);
%  for mm = 1:numel(unassigned)
%    tmp = nodes(unassigned(mm), :);
%    all6 = [[max(1,tmp(1)-1) tmp(2:3)];[min(stackSize(1),tmp(1)+1) tmp(2:3)]; [tmp(1) max(1,tmp(2)-1) tmp(3)]; [tmp(1) min(stackSize(2),tmp(2)+1) tmp(3)]];
%    all6 = [all6; [tmp(1:2) max(1,tmp(3)-1)]; [tmp(1:2) min(stackSize(3),tmp(3)+1)]];
%    tmpind = sub2ind(stackSize, all6(:,1), all6(:,2), all6(:,3));
%    tmpassign = setdiff(unique(svSeg(tmpind)), 0, 'stable');
%    if numel(tmpassign)== 1
%      groundTruth(tmpassign) = kk;
%    end
%  end
