function [nodes,edges,radii,nodeTypes,abort] = readArborTrace(fileName,validNodeTypes)
abort = false; nodes = []; edges = []; nodeTypes = [];
validNodeTypes = setdiff(validNodeTypes,1); % 1 is for soma
% read the SWC file
[nodeID, nodeType, xPos, yPos, zPos, radii, parentNodeID] = textread(fileName, '%u%d%f%f%f%f%d','commentstyle', 'shell');
% every tree should start from a node of type 1 (soma)
nodeType(find(parentNodeID==-1))=1;
% find the first soma node in the list (more than one node can be labeled as soma)
firstSomaNode = find(nodeType == 1 & parentNodeID == -1, 1);
% find the average position of all the soma nodes, and assign it as THE soma node position
somaNodes = find(nodeType == 1);
somaX = mean(xPos(somaNodes)); somaY = mean(yPos(somaNodes)); somaZ = mean(zPos(somaNodes));
somaRadius = mean(radii(somaNodes));
xPos(firstSomaNode) = somaX; yPos(firstSomaNode) = somaY; zPos(firstSomaNode) = somaZ;
radii(firstSomaNode) = somaRadius;
% change parenthood so that there is a single soma parent
parentNodeID(ismember(parentNodeID,somaNodes)) = firstSomaNode;
% delete all the soma nodes except for the firstSomaNode
nodesToDelete = setdiff(somaNodes,firstSomaNode);
nodeID(nodesToDelete)=[]; nodeType(nodesToDelete)=[];
xPos(nodesToDelete)=[]; yPos(nodesToDelete)=[]; zPos(nodesToDelete)=[];
radii(nodesToDelete)=[]; parentNodeID(nodesToDelete)=[];
% reassign node IDs due to deletions
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk)
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
% of all the nodes, retain the ones indicated in validNodeTypes
% ensure connectedness of the tree if a child is marked as valid but not some of its ancestors
validNodes = nodeID(ismember(nodeType,validNodeTypes));
additionalValidNodes = [];
for kk = 1:numel(validNodes)
  thisParentNodeID = parentNodeID(validNodes(kk)); thisParentNodeType = nodeType(thisParentNodeID);
  while ~ismember(thisParentNodeType,validNodeTypes)
    if thisParentNodeType == 1
      break;
    end
    additionalValidNodes = union(additionalValidNodes, thisParentNodeID); nodeType(thisParentNodeID) = validNodeTypes(1);
    thisParentNodeID = parentNodeID(thisParentNodeID); thisParentNodeType = nodeType(thisParentNodeID);
  end
end
% retain the valid nodes only
% the soma node is always a valid node to ensure connectedness of the tree
validNodes = [firstSomaNode; validNodes; additionalValidNodes']; validNodes = unique(validNodes);
nodeID = nodeID(validNodes); nodeType = nodeType(validNodes); parentNodeID = parentNodeID(validNodes);
xPos = xPos(validNodes); yPos = yPos(validNodes); zPos = zPos(validNodes); radii = radii(validNodes);
% reassign node IDs after deletions
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk)
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
% return the resulting tree data
nodes = [xPos yPos zPos];
edges = [nodeID parentNodeID];
edges(any(edges==-1,2),:) = [];
nodeTypes = unique(nodeType)';
