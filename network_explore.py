import json
import numpy as np
import networkx as nx
import cobra
from matplotlib import pyplot as plt


def gen_sorted_node_list(G):
	node_dict = dict(G.degree)
	to_return = sorted(node_dict,key=node_dict.get,reverse=True)
	return to_return


def reduce_graph(graphs):
	key_nodes = []
	for graph in graphs:
		sorted_nodes = gen_sorted_node_list(graph)
		accounted_for = set()
		nodes = set(graph.nodes)
		while accounted_for!=nodes:
			node = sorted_nodes.pop(0)
			key_nodes.append(node)
			neighbours = set(graph[node].keys())
			accounted_for.update(neighbours)
			accounted_for.update([node]	)
			graph.remove_node(node)
			sorted_nodes = gen_sorted_node_list(graph)
			# for neighbour in neighbours:
			# print('neighbour',neighbour)
			print('nodes', nodes)
			print('accounted_for', accounted_for)
			print('key_nodes',key_nodes)
			print('node',node)
			print('sorted_nodes',sorted_nodes)
				# sorted_nodes.remove(neighbour)
	return key_nodes

filename = 'new_data/data_hp_succ.json'

data = json.load(open(filename))

h_p = cobra.io.load_json_model(data['model'])
data = data['network']
A = np.array(data['A'])
beta1 = np.array(data['beta1']).flatten()
beta2 = np.array(data['beta2']).flatten()

print(A)

nodes = 30

print(np.shape(beta1))
print(np.shape(beta2))
print(np.shape(A))

maxes1 = np.argpartition(beta1,-nodes)[-nodes:]
small1 = A[maxes1][:,maxes1]
G = nx.from_numpy_matrix(small1)

graphs = list(nx.connected_component_subgraphs(G))

imp_nodes1 = reduce_graph(graphs)

print(imp_nodes1)

node_set = set(imp_nodes1)

node_colours = [ 'c' if x in node_set else 'r' for x in G.nodes]

plt.figure(1)
nx.draw_kamada_kawai(G,with_labels=True,font_color='k',node_color=node_colours)

maxes2 = np.argpartition(beta2,-nodes)[-nodes:]
small2 = A[maxes2][:,maxes2]
G = nx.from_numpy_matrix(small2)

graphs = list(nx.connected_component_subgraphs(G))

imp_nodes2 = reduce_graph(graphs)

print(imp_nodes2)

node_set = set(imp_nodes2)

node_colours = [ 'c' if x in node_set else 'r' for x in G.nodes]

plt.figure(2)
nx.draw_kamada_kawai(G,with_labels=True,font_color='k',node_color=node_colours)

plt.show()

# for graph in graphs:
# 	# plt.figure(i)
# 	# nx.draw_kamada_kawai(graph,with_labels=True,font_color='k')
# 	i += 1

# plt.show()