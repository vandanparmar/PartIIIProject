import network_regression
import pareto
from tqdm import tqdm
import cobra
import json
import pareto_reconstruction
from matplotlib import pyplot as plt

obj2_arr = ['RBFSa','PYRt2','PABB','O2t','NO3t2','NH4t','NAt3_1','SUCCt2r','ACt2r']

obj1_arr = ['BIOMASS_HP_published']*len(obj2_arr)

model_arr = ['h_pylori.json']*len(obj2_arr)

filename_arr = ['hp_ribo','hp_pyro','hp_fola','hp_o2_trans','hp_nit_trans','hp_amm_trans','hp_na_trans','hp_succ_neg','hp_acet_neg']

xlabel_arr = ['Riboflavin Synthesis','Pyruvate Transport','Folate Synthesis','Oxygen Transport','Nitrate Transport','Ammonia Transport','Sodium Transport','Succinate Transport','Acetate Transport']

data_str = 'new_data/data_'
figure_str = 'figures/'

lambd = 1.0
alpha = 0.0
cutoff = 95
generations = 30
individuals = 200
nodes = 20
recon_points = 200
# generations = 10
# individuals = 10
# cutoff = 70

for obj1_str,obj2_str,model_str,filename,xlabel in tqdm(list(zip(obj1_arr, obj2_arr, model_arr,filename_arr,xlabel_arr))):
	###### Pareto ##########
	network = True
	model = cobra.io.load_json_model(model_str)
	obj1 = model.reactions.get_by_id(obj1_str).flux_expression
	if "neg" in filename:
		obj2 = -model.reactions.get_by_id(obj2_str).flux_expression
		tqdm.write('Negative')
	else:
		obj2 = model.reactions.get_by_id(obj2_str).flux_expression
		tqdm.write('Positive')
	pops, vals, pareto_this = pareto.pareto(generations, individuals,model,obj1,obj2,batch = tqdm.write)

	to_save = {'obj1_str': str(obj1), 'obj2_str': str(obj2), 'model' : model_str,'pareto': [{'obj1': p.fitness.values[0], 'obj2':p.fitness.values[1],'gene_set': list(p)} for p in pareto_this]}
	pareto.plot_pareto(to_save['pareto'], 'pareto_'+filename, "H.Pylori Pareto Front", xlabel, figure_str)
	with open(data_str+filename+'.json', 'w') as outfile:
	    json.dump(to_save, outfile)
	tqdm.write("Pareto Generated.")

	##### Network regression #############
	try:
		network_regression.add_network_regression(data_str+filename+'.json', lambd, alpha, cutoff)
			
	except:
		tqdm.write('Network Regression Failed')	
		network=False
	##### Pareto Reconstruction ########
	if network==True:
		pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x = pareto_reconstruction.reconstruct(data_str+filename+'json', nodes, recon_points, model, obj1, obj2)
		plt.clf()
		plt.plot(pareto_x,pareto_y,'*',color='k',label='Pareto Optimal')
		plt.plot(pareto_left[:,1], pareto_left[:,0],'r.',label='Left')
		plt.plot(pareto_right[:,1], pareto_right[:,0],'c.',label='Right')
		plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.',label='Noise')
		plt.xlabel(xlabel)
		plt.ylabel('Biomass Production')
		plt.title('H. Pylori Pareto Front Reconstruction')
		plt.legend()
		plt.savefig('recon_'+figure_str+filename+'.png')
		plt.savefig('recon_'+figure_str+filename+'.eps')
		plt.clf()

