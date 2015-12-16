function subgraphNodes = randomlyPartitionBinaryGraph(binsaff, maxSubgraphSize)

availableNodes                 = 1:size(binsaff, 1);
subgraphNodes                  = cell(0);
subgraphSizes                  = [];
while true
  thisSeed                     = randi(numel(availableNodes));
  tmp                          = growSubgraphOnBinaryGraph(binsaff(availableNodes, availableNodes), thisSeed, maxSubgraphSize);
  subgraphNodes{end+1}         = availableNodes(tmp);
  subgraphSizes(end+1)         = numel(tmp);
  availableNodes               = availableNodes(setdiff(1:numel(availableNodes), tmp, 'stable'));
  if numel(availableNodes)<maxSubgraphSize
    subgraphNodes{end+1}       = availableNodes;
    subgraphSizes(end+1)       = numel(availableNodes);
    break;
  end
end

[sizes, sorted]                = sort(subgraphSizes);
while sizes(1)<maxSubgraphSize & numel(sizes)>1
  pos                          = find(cumsum(sizes)>maxSubgraphSize, 1);             
  if isempty(pos)
    pos                        = numel(sizes);
  end
  subgraphNodes{sorted(1)}     = cell2mat(subgraphNodes(sorted(1:pos)));
  subgraphSizes(sorted(1))     = numel(subgraphNodes{sorted(1)});
  subgraphNodes(sorted(2:pos)) = [];
  subgraphSizes(sorted(2:pos)) = [];
  [sizes, sorted]              = sort(subgraphSizes);
end
disp(cellfun(@numel, subgraphNodes));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeSubset = growSubgraphOnBinaryGraph(binsaff, node, maxSubgraphSize)

nodeSubset      = node;
newlyAdded      = node;
oldCount        = 1;
while true
  oldNodeSubset = nodeSubset;
  nodeSubset    = union(nodeSubset, find(any(binsaff(newlyAdded, :), 1)), 'stable');
  newlyAdded    = setdiff(nodeSubset, oldNodeSubset, 'stable');
  count         = numel(nodeSubset);
  if count==oldCount | count>=maxSubgraphSize
    break;
  else
    oldCount    = count;
  end
end
