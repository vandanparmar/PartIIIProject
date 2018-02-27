import json
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import cobra


data = json.load(open('networks2.json'))

h_p = cobra.io.load_json_model('iJO1366.json')

A = np.array(data['A'])
beta1 = np.array(data['beta1']).flatten()
beta2 = np.array(data['beta2']).flatten()

nodes = 30

print(np.shape(beta1))
print(np.shape(beta2))
print(np.shape(A))

plt.figure(1)
maxes1 = np.argpartition(beta1,-nodes)[-nodes:]
small1 = A[maxes1][:,maxes1]
G = nx.from_numpy_matrix(small1)

gene_list = list(map(lambda i : h_p.genes[i],maxes1))

mapping=dict(zip(G.nodes(),gene_list))
G = nx.relabel_nodes(G,mapping)
nx.draw_kamada_kawai(G,with_labels=True,font_color='k')
for i in gene_list:
	print(i,i.reactions)


plt.figure(2)
maxes2 = np.argpartition(beta2,-nodes)[-nodes:]
small2 = A[maxes2][:,maxes2]
G = nx.from_numpy_matrix(small2)
mapping=dict(zip(G.nodes(),maxes2))
G = nx.relabel_nodes(G,mapping)
nx.draw_kamada_kawai(G,with_labels=True,font_color='k')
plt.show()
