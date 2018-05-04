import json
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import cobra
import network_regression as nr

data = json.load(open('saturated_data/ratios/datahp_no_loop_o2.json'))

filename1 = 'saturated_data/ratios/datahp_no_loop_o2.json'
filename2 = 'saturated_data/ratios/datahp_no_loop_co2.json'
lambd = 0.1
alpha = 0.1
cutoff = 99.99
lambdas = [0.0,0.01,0.1,1]
h_p = cobra.io.load_json_model('h_pylori.json')
# for i,lambd in enumerate(lambdas):
	# betal,betar = nr.add_network_regression(filename1, lambd, alpha, cutoff)
betal,betar,W,A,S = nr.add_network_regression(filename2, lambd, alpha, cutoff)
plt.imshow(S)	
plt.show()
# 	plt.figure(i+1)
# 	plt.hist(betar.flatten().T,histtype='step',bins = 100,label='right'+str(lambd),density = False,log = False)
# 	plt.hist(betal.flatten().T,histtype='step',bins = 100,label ='left'+str(lambd),density = False,log = False)
# plt.legend()
# plt.show()

# A = np.array(data['A'])
# beta1 = np.array(data['beta1']).flatten()
# beta2 = np.array(data['beta2']).flatten()

# nodes = 30

# print(np.shape(beta1))
# print(np.shape(beta2))
# print(np.shape(A))

# plt.figure(1)
# maxes1 = np.argpartition(beta1,-nodes)[-nodes:]
# small1 = A[maxes1][:,maxes1]
# G = nx.from_numpy_matrix(small1)

# gene_list = list(map(lambda i : h_p.genes[i],maxes1))

# mapping=dict(zip(G.nodes(),gene_list))
# G = nx.relabel_nodes(G,mapping)
# nx.draw_kamada_kawai(G,with_labels=True,font_color='k')
# for i in gene_list:
# 	print(i,i.reactions)


# plt.figure(2)
# maxes2 = np.argpartition(beta2,-nodes)[-nodes:]
# small2 = A[maxes2][:,maxes2]
# G = nx.from_numpy_matrix(small2)
# mapping=dict(zip(G.nodes(),maxes2))
# G = nx.relabel_nodes(G,mapping)
# nx.draw_kamada_kawai(G,with_labels=True,font_color='k')
# plt.show()
