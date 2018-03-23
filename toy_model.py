import cobra
import pareto
import numpy as np
from matplotlib import pyplot as plt
import kink_finder
from tqdm import tqdm

model = cobra.Model('model')



a = cobra.Metabolite('a')
b = cobra.Metabolite('b')
c = cobra.Metabolite('c')
d = cobra.Metabolite('d')
e = cobra.Metabolite('e')
f = cobra.Metabolite('f')
g = cobra.Metabolite('g')
h = cobra.Metabolite('h')

#metabolites = ['a','b','c','d']
metabolites = ['b','c']


#intake_genes = ['','','','']
intake_genes = ['g1','g2']


# reaction_names = ['biomass','s1','s2','s3','t1','t2']
# reaction_rules = ['g1','g2','g2','g1','(g1 or g2)','(g1 or g2)']

reaction_names = ['biomass','s1','s2','s3','t1']
reaction_rules = ['g1','g2','g3','g2','g4']

stochiometries = [{a:1.0,c:-1.0},{d:1.0,c:-1.0},{a:-1.0,b:1.0},{b:1.0,d:-1.0},{a:1.0,d:-1.0}]

for metabolite,gene in zip(metabolites,intake_genes):
	reaction_i = cobra.Reaction(metabolite+'_intake')
	reaction_i.lower_bound = -3.0
	reaction_i.upper_bound = 5.0
	reaction_i.add_metabolites({
		eval(metabolite) : -1.0
		})
	reaction_i.gene_reaction_rule = gene
	model.add_reactions([reaction_i])

for reaction_name,reaction_rule,stochiometry in zip(reaction_names, reaction_rules,stochiometries):
	reaction = cobra.Reaction(reaction_name)
	if reaction_name == 'biomass':
		reaction.bounds = [-0.0,10.0]
	elif reaction_name == 't1':
		reaction.bounds = [-10.0,10.0]
	elif reaction_name == 's3':
		reaction.bounds = [5.0,8.6]
	else:
		reaction.bounds = [-100.0,100.0]
	reaction.add_metabolites(stochiometry)
	reaction.gene_reaction_rule = reaction_rule
	model.add_reactions([reaction])


generations = 30
pop_size = 1000

obj1_str = 'biomass'
obj2_str = 's3'


bio = model.reactions.get_by_id(obj1_str).flux_expression
s3 = model.reactions.get_by_id(obj2_str).flux_expression

obj1 = bio
obj2 = s3

bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))
lb = bounds[:,0]
lb = np.clip(lb, a_min=-100.0,a_max = 0.0)
ub = bounds[:,1]
ub = np.clip(ub,a_min= 0.0,  a_max=100.0)
gene_dict = {x.id:i for i,x in enumerate(model.genes)}
condts = pareto.gene_condts(model,gene_dict)


pops,vals,pareto_this =	pareto.pareto(generations,pop_size,model, obj1, obj2)


points = list(zip(*[p.fitness.values for p in pareto_this]))

c1 = np.array([194, 45, 213])/255
c2 = np.array([52, 188, 110])/255

def c_map(c1,c2,v_min,v_max):
	diff1 = c1-c2
	diff2 = v_max-v_min
	ratio = diff1/diff2
	return (lambda i: tuple(ratio*(i-v_min)+c2))


fig = plt.figure()
ax = fig.add_subplot(111)
tot = len(vals)
c_map_i = c_map(c1,c2,0,tot)
print(tot)
for i,val in enumerate(vals):

	ax.scatter(val[1],val[0],marker='*',c=c_map_i(i))
	ax.plot(val[1],val[0],'r--',linewidth=0.5,c=c_map_i(i))
plt.ylabel('$f_1$')
plt.xlabel('$f_2$')
plt.title('Pareto')
plt.show()





to_save = {'pareto': [{'obj1': p.fitness.values[0], 'obj2':p.fitness.values[1],'gene_set': list(p)} for p in pareto_this]}
pareto_this = to_save['pareto']
pareto_this = sorted(pareto_this,key= lambda x : x['obj2'])
print(pareto_this[0]['obj1'])
print(pareto_this[-1]['obj1'])

y_plot = np.array(list(map(lambda i : i['obj1'],pareto_this)))
x_plot = np.array(list(map(lambda i : i['obj2'],pareto_this)))
genes = np.array(list(map(lambda i : i ['gene_set'],pareto_this)))

# const = np.arange(0, len(model.reactions)*0.05, 0.05)
x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot,y_plot)
plt.subplot(311)
plt.xlabel('s3 Flux')
plt.ylabel('Biomass Production')
plt.plot(x_plot,y_plot,c='b', label = 'Pareto Points',marker='.',linestyle='None')
plt.plot(x0,y0,'*',c='r', label = 'Phase Transition Point')
plt.plot(x_plot,kink_finder.piecewise_linear(x_plot,x0,y0,k1,k2),c='lawngreen', label = 'Fitted Line')
plt.legend()


fluxes = list(map(lambda x: pareto.long_evaluate(x['gene_set'],obj1,obj2,model,condts,lb,ub).values,tqdm(pareto_this)))
# fluxes += const
reaction_names = list(map(lambda x: x.id, model.reactions))
plt.subplot(312)
fluxes = np.array(fluxes)
line_objs = plt.plot(x_plot,fluxes)
plt.legend(line_objs,reaction_names)

plt.subplot(313)
line_objs = plt.plot(x_plot,genes,'.',markersize=2)
gene_names = np.array(list(map(lambda x : x.id,model.genes)))
plt.legend(line_objs,gene_names)
plt.show()
