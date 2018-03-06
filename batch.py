import network_regression
import pareto
from tqdm import tqdm
import cobra
import json

# obj2_arr = ['SUCCt2_2pp','ACt4pp','EX_o2_e','SUCCt2r','ACt2r','EX_o2_e']



# obj1_arr = ['BIOMASS_Ec_iJO1366_WT_53p95M']*3
# obj1_arr += ['BIOMASS_HP_published']*3

# model_arr = ['iJO1366.json']*3
# model_arr += ['h_pylori.json']*3

# filename_arr = ['ecoli_succ','ecoli_acet','ecoli_o2','hp_succ','hp_acet','hp_o2']


obj2_arr = ['ACt2r']



# obj1_arr = ['BIOMASS_Ec_iJO1366_WT_53p95M']*3
obj1_arr = ['BIOMASS_HP_published']

# model_arr = ['iJO1366.json']*
model_arr = ['h_pylori.json']

filename_arr = ['hp_acet']


data_str = 'new_data/data_'
figure_location = 'figures/'

lambd = 1.0
alpha = 0.0
cutoff = 95
generations = 30
individuals = 200


# generations = 10
# individuals = 10
# cutoff = 70


for obj1_str,obj2_str,model_str,filename in tqdm(list(zip(obj1_arr, obj2_arr, model_arr,filename_arr))):
	###### Pareto ##########

	model = cobra.io.load_json_model(model_str)
	obj1 = model.reactions.get_by_id(obj1_str).flux_expression
	obj2 = model.reactions.get_by_id(obj2_str).flux_expression

	pops, vals, pareto_this = pareto.pareto(generations, individuals,model,obj1,obj2,batch = tqdm.write)

	to_save = {'obj1_str': str(obj1), 'obj2_str': str(obj2), 'model' : model_str,'pareto': [{'obj1': p.fitness.values[0], 'obj2':p.fitness.values[1],'gene_set': list(p)} for p in pareto_this]}

	with open(data_str+filename+'.json', 'w') as outfile:
	    json.dump(to_save, outfile)
	tqdm.write("Pareto Generated.")

	##### Network regression #############
	try:
		network_regression.add_network_regression(data_str+filename+'.json', lambd, alpha, cutoff)
	except:
		tqdm.write('Network Regression Failed')	