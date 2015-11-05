from snapvx import *
from numpy import genfromtxt
from joblib import Parallel, delayed
import numpy as np
import csv
import time
#Helper function for node objective
#Takes in a row from the CSV file, returns an optimization problem
# def node_obj(data):
#   x = Variable(1,name='x');
#   return power(norm(x - float(data[0])),2);

# #Helper function for edge objective
# def laplace_reg(src, dst, data):
#   weight = weight_dict[(src,dst)];
#   return weight*norm(src['x'] - dst['x']);

# gvx = TGraphVX()
# my_data = genfromtxt('superVoxelMeans.csv', delimiter=',')
# for i in range(my_data.shape[0]):
#   a = my_data[i,];
#   x = Variable(my_data.shape[1],name='x');
#   obj = square(norm(x-a));
#   gvx.AddNode(i+1, Objective=obj);




my_data = genfromtxt('edges_weights.csv', delimiter=',')
weight_dict = {};
unique_edges = set([]);
for i in range(my_data.shape[0]):
  srcID = int(my_data[i,0]);  desID = int(my_data[i,1]);
  elements = frozenset([srcID,desID]);
  unique_edges.add(elements);
  weight = my_data[i,2];
  weight_dict[(srcID,desID)] = weight;

unique_edges_list = list(unique_edges)

#Load in Edge List to build graph with default node/edge objectives
gvx = LoadEdgeList('edges_supervoxel.txt');

my_data = genfromtxt('superVoxelMeans.csv', delimiter=',');
SVsize = genfromtxt('superVoxelSizes.csv', delimiter=',');
for node in gvx.Nodes():
  a = my_data[node.GetId()-1,];
  x = Variable(my_data.shape[1],name='x');
  gvx.SetNodeObjective(node.GetId(), square(norm(x-a)));#SVsize[node.GetId()-1]*


#Bulk Load node objectives:
#Takes one row of the CSV, uses that as input to node_obj
#There is also an (optional) input of specifying which nodes each row of the CSV refers to

#Bulk load edge objectives for all edges

# def add_edge_objective(edges_list):
#   count = 0;
#   for item in edges_list:
#     edge = set(item); 
#     src = edge.pop(); dst = edge.pop();
#     print 'source = '+repr(src)+' des = '+repr(dst)
#     if (src,dst) in weight_dict.keys():
#       weight = weight_dict[(src,dst)];
#     else:
#       weight = weight_dict[(dst,src)];
#     obj = norm(gvx.GetNodeVariables(src)['x'] - gvx.GetNodeVariables(dst)['x'])
#     gvx.SetEdgeObjective(src, dst, weight*obj);
#     count = count+1;
#     if not (count%100):
#       print 'Time taken = '+repr(time.time()-start_time)
#       start_time = time.time()




# num_cores = 8;
# each_core_count = int(len(unique_edges_list)/num_cores);
# core_index = [];
# for i in range(num_cores):
#   if (i==(num_cores-1)):
#     core_index.append(np.arange(i*each_core_count, len(unique_edges_list)))
#   else:
#     core_index.append(np.arange(i*each_core_count, (i+1)*each_core_count))



# result = (Parallel(n_jobs=num_cores)(delayed(add_edge_objective)([unique_edges_list[a] for a in core_index[i]]) for i in range(num_cores)))
# notusefulIDpart = [];
# for temp in result:
#   notusefulIDpart.extend(temp)
#lambda_list = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2]

#for lambda_value in lambda_list:
lambda_value = 1;
print 'Solving for lambda = '+repr(lambda_value)
count = 0;
start_time = time.time();
for edge in gvx.Edges():
  src = edge.GetSrcNId();
  dst = edge.GetDstNId();
  if (src,dst) in weight_dict.keys():
    weight = weight_dict[(src,dst)];#lambda_value*
  else:
    weight = weight_dict[(dst,src)];#lambda_value*
  obj = norm(gvx.GetNodeVariables(src)['x'] - gvx.GetNodeVariables(dst)['x'])
  gvx.SetEdgeObjective(src, dst, weight*obj);
  count = count+1;
  if not (count%100):
    print 'Time taken = '+repr(time.time()-start_time)
    start_time = time.time()





# for key in weight_dict.keys():
#     weight = weight_dict[key];
#     src = key[0]; dst = key[1];
#     #x = Variable(my_data.shape[1],name='x');
#     obj = norm(src['x'] - dst['x']);
#     gvx.AddEdge(key[0], key[1], Objective = weight*obj);


gvx.Solve(UseADMM=True, Verbose=True);
with open('nodeValue_'+str(lambda_value)+'.csv', 'w', newline='') as csvfile:
  rewriter = csv.writer(csvfile, delimiter=',',quotechar='''\"''', quoting=csv.QUOTE_MINIMAL)
  for node in gvx.Nodes():
    rewriter.writerow(np.append([node.GetId()],gvx.GetNodeValue(node.GetId(),'x')));
