import numpy as np
import json
import pareto
import cobra
from matplotlib import pyplot as plt
from tqdm import tqdm


nodes = 30
recon_points = 100
sp_max = 2
sp_min = 0


# model_str = data['model']
# network = data['network']
# pareto_this = data['pareto']


def eval_pareto(model ,obj1 ,obj2 ,gene_set):
	evaluate = pareto.evaluate
	bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))
	lb = bounds[:,0]
	lb = np.clip(lb, a_min=-100.0,a_max = 0.0)
	ub = bounds[:,1]
	ub = np.clip(ub,a_min= 0.0,  a_max=100.0)
	gene_dict = {x.id:i for i, x in enumerate(model.genes)}
	condts = pareto.gene_condts(model, gene_dict)
	# evaluate(individual,obj1,obj2,model,condts,lb,ub)
	gene_set = tqdm(gene_set)
	to_return = list(map(lambda i : evaluate(i, obj1, obj2, model, condts, lb, ub), gene_set))
	return np.array(to_return)


network = json.load(open('networks_ecoli_succ.json'))
model_str = 'iJO1366.json'
pareto_data = json.load(open('data_big_succinate.json'))
obj1_str = 'BIOMASS_Ec_iJO1366_WT_53p95M'
obj2_str = 'SUCCt2_2pp'

y_plot = np.array(list(map(lambda i : i['obj1'], pareto_data)))
x_plot = np.array(list(map(lambda i : i['obj2'], pareto_data)))
pareto_genes = np.array(list(map(lambda i : i['gene_set'], pareto_data)))




h_p = cobra.io.load_json_model(model_str)

obj1 = h_p.reactions.get_by_id(obj1_str).flux_expression
obj2 = h_p.reactions.get_by_id(obj2_str).flux_expression


print(obj1)
print(obj2)

beta1 = np.array(network['beta1']).flatten()
beta2 = np.array(network['beta2']).flatten()


maxes1 = np.argpartition(beta1, -nodes)[-nodes:]
maxes2 = np.argpartition(beta2, -nodes)[-nodes:]

print(np.shape(maxes1))
print(np.shape(maxes2))

n_genes = len(h_p.genes)
to_pareto = np.random.uniform(low=0.0, high=2.0, size=(recon_points*2, n_genes))
to_pareto_noise = to_pareto
# to_set = np.linspace(sp_min, sp_max, recon_points)
# to_set = np.reshape(np.repeat(to_set, nodes), (-1, nodes))

to_set = np.ones((recon_points,nodes))
to_set -= 0.7

to_pareto[::2, maxes1] = to_set
to_pareto[1::2, maxes2] = to_set

print(to_pareto)

# plt.hist(pareto_genes.flatten(),histtype='step')
# plt.hist(to_pareto.flatten(),histtype='step')
# plt.show()

pareto_new = eval_pareto(h_p, obj1, obj2, to_pareto)
pareto_noise = eval_pareto(h_p,obj1,obj2,to_pareto_noise)
print(pareto_new)

print(np.shape(pareto_new))



plt.plot(pareto_new[1::2,1], pareto_new[1::2,0],'c.')
plt.plot(pareto_new[::2,1], pareto_new[::2,0],'r.')
plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.')
# plt.plot(x_plot,y_plot,'r.')
plt.show()

