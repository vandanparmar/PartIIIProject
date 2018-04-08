import numpy as np
import json
import pareto
import cobra
from matplotlib import pyplot as plt
from tqdm import tqdm
from itertools import repeat
import networkx as nx
import kink_finder
import network_explore
import multiprocessing

# model_str = data['model']
# network = data['network']
# pareto_this = data['pareto']

def show_network(maxes,A):
	small = A[maxes][:,maxes]
	G = nx.from_numpy_matrix(small)
	mapping=dict(zip(G.nodes(),maxes))
	G = nx.relabel_nodes(G,mapping)
	nx.draw_spring(G,with_labels=True,font_color='k')
	plt.show()
	return

def prep_fba_set(points, descriptor):
	to_return = list(zip(points[:,0].tolist(),points[:,1].tolist(),repeat(descriptor)))
	return to_return


def find_optimal(fba_points):
	fba_points.sort(key = lambda x : x[0],reverse = True)
	to_return = []
	current_max = -np.inf
	for point in fba_points:
		if point[1] > current_max:
			current_max = point[1]
			to_return.append(point)
	return to_return


def eval_pareto(model ,obj1 ,obj2 ,gene_set,cores):
	bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))
	lb = bounds[:,0]
	lb = np.clip(lb, a_min=-np.inf,a_max = np.inf)
	ub = bounds[:,1]
	ub = np.clip(ub,a_min= -np.inf,  a_max=np.inf)
	gene_dict = {x.id:i for i, x in enumerate(model.genes)}
	condts = pareto.gene_condts(model, gene_dict)
	# evaluate(individual,obj1,obj2,model,condts,lb,ub)
	if cores==0:
		evaluate = pareto.evaluate
		gene_set = tqdm(gene_set)
		to_return = list(map(lambda i : evaluate(i, obj1, obj2, model, condts, lb, ub), gene_set))
	else:
		evaluate = pareto.eval_wrapper
		pool = multiprocessing.Pool(cores)	
		args = zip(gene_set, repeat(obj1),repeat(obj2),repeat(model),repeat(condts),repeat(lb),repeat(ub))
		to_return = list(tqdm(pool.imap(evaluate, args,chunksize=50),total = len(gene_set)))
	return np.array(to_return)

def origin_reconstruct(filename,h_p,obj1,obj2):
	data = json.load(open(filename))
	pareto_data = data['pareto']
	pareto_genes = np.array(list(map(lambda i : i['gene_set'], pareto_data)))[::4]
	pareto = eval_pareto(h_p, obj1, obj2, pareto_genes)
	return pareto

def create_network_vals(key_nodes,G,maxes,recon_points,low,high,pareto_genes):
	node_dict = {v:i for i,v in enumerate(maxes)}
	key_node_dict = {v:i for i,v in enumerate(key_nodes)}	
	small_indices = np.random.random_integers(low, high, size=(recon_points,len(key_nodes)))
	indices = np.zeros((recon_points,len(maxes)))	
	for i,maxi in enumerate(maxes):
		if node_dict[maxi] in key_nodes:
			indices[:,i] = small_indices[:,key_node_dict[node_dict[maxi]]]
		else:
			neighbours = list(G[node_dict[maxi]].keys())
			neighbours = [x for x in neighbours if x in key_nodes]
			neighbours = list(map(lambda i : key_node_dict[i],neighbours))
			choices = np.random.choice(neighbours, recon_points)
			indices[:,i] = small_indices[np.arange(recon_points),choices]
	indices = np.array(indices, dtype=np.int)
	to_set = pareto_genes[indices,maxes]
	return to_set

def network_reconstruct(filename,nodes,recon_points,h_p,obj1,obj2,cores=0):
	data = json.load(open(filename))
	network = data['network']
	pareto_data = data['pareto']
	y_plot = np.array(list(map(lambda i : i['obj1'], pareto_data)))
	x_plot = np.array(list(map(lambda i : i['obj2'], pareto_data)))
	pareto_genes = np.array(list(map(lambda i : i['gene_set'], pareto_data)))

	bounds = np.array(list(map(lambda reaction : reaction.bounds, h_p.reactions)))

	x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot, y_plot)
	phase_trans = np.abs(x_plot-x0).argmin()

	beta1 = np.array(network['beta1']).flatten()
	beta2 = np.array(network['beta2']).flatten()

	maxes1 = np.argpartition(beta1, -nodes)[-nodes:]
	maxes2 = np.argpartition(beta2, -nodes)[-nodes:]

	n_genes = len(h_p.genes)
	to_pareto = np.random.uniform(low=0.0, high=2.0, size=(recon_points*2, n_genes))
	to_pareto_noise = np.random.uniform(low=0.0, high=2.0, size=(recon_points,n_genes))

	A = np.array(network['A'])

	small1 = A[maxes1][:,maxes1]
	small2 = A[maxes2][:,maxes2]

	G1 = nx.from_numpy_matrix(small1)
	G2 = nx.from_numpy_matrix(small2)
	
	graphs1 = list(nx.connected_component_subgraphs(G1))
	graphs2 = list(nx.connected_component_subgraphs(G2))

	key_nodes1 = network_explore.reduce_graph(graphs1)
	key_nodes2 = network_explore.reduce_graph(graphs2)

	node_colours = [ 'c' if x in key_nodes1 else 'r' for x in G1.nodes]

	plt.figure(3)
	nx.draw_kamada_kawai(G1,with_labels=True,font_color='k',node_color=node_colours)

	# plt.show()

	G1 = nx.from_numpy_matrix(small1)
	G2 = nx.from_numpy_matrix(small2)

	to_set_1 = create_network_vals(key_nodes1, G1, maxes1, recon_points, low=0, high=phase_trans, pareto_genes=pareto_genes)
	to_set_2 = create_network_vals(key_nodes2, G2, maxes2, recon_points, low=phase_trans+1,high=len(x_plot)-1, pareto_genes=pareto_genes)
	to_pareto[::2, maxes1] = to_set_1
	to_pareto[1::2, maxes2] = to_set_2

	pareto_new = eval_pareto(h_p, obj1, obj2, to_pareto,cores)

	for i,x in enumerate(bounds):
		h_p.reactions[i].bounds = x

	pareto_noise = eval_pareto(h_p,obj1,obj2,to_pareto_noise,cores)

	for i,x in enumerate(bounds):
		h_p.reactions[i].bounds = x

	pareto_left = pareto_new[::2]
	pareto_right = pareto_new[1::2]

	paretos = prep_fba_set(pareto_left, 'left')
	paretos.extend(prep_fba_set(pareto_right,'right'))
	paretos.extend(prep_fba_set(pareto_noise,'noise'))

	pareto_collection = find_optimal(paretos)
	pareto_y = list(map(lambda x : x[0],pareto_collection))
	pareto_x = list(map(lambda x : x[1],pareto_collection))

	to_add = {'pareto_left':pareto_left.tolist(),
		'pareto_right':pareto_right.tolist(),
		'pareto_noise':pareto_noise.tolist(),
		'pareto_y':pareto_y,
		'pareto_x':pareto_x }
	data['network_recon'] = to_add
	with open(filename, 'w') as outfile:
	    json.dump(data, outfile)

	return pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x



def reconstruct(filename,nodes,recon_points,h_p,obj1,obj2,cores=0):
	data = json.load(open(filename))
	network = data['network']
	pareto_data = data['pareto']
	y_plot = np.array(list(map(lambda i : i['obj1'], pareto_data)))
	x_plot = np.array(list(map(lambda i : i['obj2'], pareto_data)))
	pareto_genes = np.array(list(map(lambda i : i['gene_set'], pareto_data)))

	bounds = np.array(list(map(lambda reaction : reaction.bounds, h_p.reactions)))

	x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot, y_plot)
	phase_trans = np.abs(x_plot-x0).argmin()

	beta1 = np.array(network['beta1']).flatten()
	beta2 = np.array(network['beta2']).flatten()

	maxes1 = np.argpartition(beta1, -nodes)[-nodes:]
	maxes2 = np.argpartition(beta2, -nodes)[-nodes:]

	n_genes = len(h_p.genes)
	to_pareto = np.random.uniform(low=0.0, high=2.0, size=(recon_points*2, n_genes))
	to_pareto_noise = np.random.uniform(low=0.0, high=2.0, size=(recon_points,n_genes))

	indices1 = np.random.random_integers(low=0,high=phase_trans,size=(recon_points,nodes))
	indices2 = np.random.random_integers(low=phase_trans+1,high=len(x_plot)-1,size=(recon_points,nodes))

	to_set_1 = pareto_genes[indices1,maxes1]
	to_set_2 = pareto_genes[indices2,maxes2]
	to_pareto[::2, maxes1] = to_set_1
	to_pareto[1::2, maxes2] = to_set_2

	pareto_new = eval_pareto(h_p, obj1, obj2, to_pareto,cores)

	for i,x in enumerate(bounds):
		h_p.reactions[i].bounds = x

	pareto_noise = eval_pareto(h_p,obj1,obj2,to_pareto_noise,cores)

	pareto_left = pareto_new[::2]
	pareto_right = pareto_new[1::2]

	paretos = prep_fba_set(pareto_left, 'left')
	paretos.extend(prep_fba_set(pareto_right,'right'))
	paretos.extend(prep_fba_set(pareto_noise,'noise'))

	pareto_collection = find_optimal(paretos)
	pareto_y = list(map(lambda x : x[0],pareto_collection))
	pareto_x = list(map(lambda x : x[1],pareto_collection))

	to_add = {'pareto_left':pareto_left.tolist(),
		'pareto_right':pareto_right.tolist(),
		'pareto_noise':pareto_noise.tolist(),
		'pareto_y':pareto_y,
		'pareto_x':pareto_x }
	data['recon'] = to_add
	with open(filename, 'w') as outfile:
	    json.dump(data, outfile)

	return pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x

if __name__ == '__main__':
	nodes = 40
	recon_points = 300

	data = json.load(open('new_data/data_hp_acet.json'))
	network = data['network']
	model_str = data['model']
	pareto_data = data['pareto']
	obj1_str = data['obj1_str']
	obj2_str = data['obj2_str']

	# obj1_str = 'BIOMASS_Ec_iJO1366_WT_53p95M'

	obj1_str = 'BIOMASS_HP_published'
	obj2_str = 'ACt2r'
	# obj2_str = 'EX_o2_e'

	y_plot = np.array(list(map(lambda i : i['obj1'], pareto_data)))
	x_plot = np.array(list(map(lambda i : i['obj2'], pareto_data)))
	pareto_genes = np.array(list(map(lambda i : i['gene_set'], pareto_data)))

	x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot, y_plot)
	phase_trans = np.abs(x_plot-x0).argmin()


	h_p = cobra.io.load_json_model(model_str)

	obj1 = h_p.reactions.get_by_id(obj1_str).flux_expression
	obj2 = h_p.reactions.get_by_id(obj2_str).flux_expression


	print(obj1)
	print(obj2)

	beta1 = np.array(network['beta1']).flatten()
	beta2 = np.array(network['beta2']).flatten()


	maxes1 = np.argpartition(beta1, -nodes)[:nodes]
	maxes2 = np.argpartition(beta2, -nodes)[:nodes]

	print(set(maxes1)-set(maxes2))
	print(maxes1)
	print(maxes2)


	print(np.shape(maxes1))
	print(np.shape(maxes2))

	n_genes = len(h_p.genes)
	to_pareto = np.random.uniform(low=0.0, high=2.0, size=(recon_points*2, n_genes))
	to_pareto_noise = to_pareto
	# to_set = np.linspace(sp_min, sp_max, recon_points)
	# to_set = np.reshape(np.repeat(to_set, nodes), (-1, nodes))

	# to_set = np.ones((recon_points,nodes))
	# to_set += 1.0

	# plt.plot(x_plot[:phase_trans],y_plot[:phase_trans],'r.')
	# plt.plot(x_plot[phase_trans:],y_plot[phase_trans:],'g.')
	# plt.show()

	len1 = phase_trans
	len2 = len(x_plot)-phase_trans

	indices1 = np.random.random_integers(low=0,high=phase_trans,size=(recon_points,nodes))
	indices2 = np.random.random_integers(low=phase_trans+1,high=len(x_plot)-1,size=(recon_points,nodes))

	to_set_1 = pareto_genes[indices1,maxes1]
	to_set_2 = pareto_genes[indices2,maxes2]

	means1 = np.mean(pareto_genes[:phase_trans,maxes1],axis=0)
	means2 = np.mean(pareto_genes[phase_trans:,maxes2],axis=0)

	# means1 = np.mean(pareto_genes[:,maxes1],axis=0)
	# means2 = np.mean(pareto_genes[:,maxes2],axis=0)

	show_network(maxes1, np.array(network['A']))
	show_network(maxes2, np.array(network['A']))

	# to_set_1 = np.random.normal(0.0, 0.125, (recon_points,nodes))
	# to_set_1 = np.clip(to_set_1, -0.25, 0.25)
	# to_set_1 += means1
	# to_set_1 = np.clip(to_set_1, 0.0, np.inf)


	# to_set_2 = np.random.normal(0.0, 0.125, (recon_points,nodes))
	# to_set_2 = np.clip(to_set_2, -0.25, 0.25)
	# to_set_2 += means2
	# to_set_2 = np.clip(to_set_2, 0.0, np.inf)

	# plt.hist(to_set_1.flatten(),histtype='step')
	# plt.hist(to_set_2.flatten(),histtype='step')
	# plt.show()

	to_pareto[::2, maxes1] = to_set_1
	to_pareto[1::2, maxes2] = to_set_2


	for max1 in maxes2:
		plt.subplot()
		plt.hist(pareto_genes[phase_trans:,max1].flatten(),normed=True, histtype='step',bins=100)
		plt.hist(to_pareto[1::2,max1].flatten(),normed=True,linestyle='dashed',histtype='step',bins=100)
	plt.show()

	# plt.hist(to_pareto.flatten(),histtype='step')
	# plt.show()


	pareto_new = eval_pareto(h_p, obj1, obj2, to_pareto)
	pareto_noise = eval_pareto(h_p,obj1,obj2,to_pareto_noise)


	pareto_left = pareto_new[::2]
	pareto_right = pareto_new[1::2]

	paretos = prep_fba_set(pareto_left, 'left')
	paretos.extend(prep_fba_set(pareto_right,'right'))
	paretos.extend(prep_fba_set(pareto_noise,'noise'))

	pareto_collection = find_optimal(paretos)
	pareto_y = list(map(lambda x : x[0],pareto_collection))
	pareto_x = list(map(lambda x : x[1],pareto_collection))


	plt.plot(pareto_x,pareto_y,'*',color='k')

	plt.plot(pareto_new[::2,1], pareto_new[::2,0],'r.')
	plt.plot(pareto_new[1::2,1], pareto_new[1::2,0],'c.')
	plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.')
	# plt.plot(x_plot,y_plot,'r.')
	plt.show()