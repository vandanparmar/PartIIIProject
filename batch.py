import matplotlib
matplotlib.use('agg')
import network_regression
import pareto
from tqdm import tqdm
import cobra
import json
import pareto_reconstruction
from matplotlib import pyplot as plt
import numpy as np

# obj2_arr = ['EX_asp__L_e','EX_cit_e','ACt2r','G6PDH2r','GALT','G6PDH1','HBZOPT','H2Ot','H2CO3TP','GALKr','HCO3E']
# obj1_arr = ['EX_o2_e','ALAD_L','EX_o2_e','EX_arg__L_e','EX_arg__L_e','EX_arg__L_e','EX_o2_e','EX_o2_e','EX_o2_e','EX_arg__L_e','EX_o2_e']
# filename_arr = ['hp_asp_o2','hp_cit_o2','hp_acet_o2','hp_pent_arg2','hp_gala_arg','hp_pent_arg','hp_ubi_o2','hp_h2o_o2','hp_acid_o2','hp_galak_arg','hp_acide_o2']
# xlabel_arr = ['L-Aspartate Uptake','Citrate Uptake','Acetate Transport','Pentose Phosphate Transport','Galactose Synthesis','Pentose Phosphate Transport','Ubiquinone Synthesis','Water Transport','Acid Transport','Galactokinase Synthesis','Acid Equilibriation']
# ylabel_arr = ['Oxygen Uptake','Oxygen Uptake','Oxygen Uptake','Arginine Uptake','Arginine Uptake','Arginine Uptake','Oxygen Uptake','Oxygen Uptake','Oxygen Uptake','Arginine Uptake','Oxygen Uptake']
obj2_arr = []
model = cobra.io.load_json_model('h_pylori.json')
for reaction in model.reactions:
	if 'EX_' in reaction.id:
		obj2_arr.append(reaction.id)

obj2_arr = obj2_arr
obj1_arr = ['BIOMASS_HP_published']*len(obj2_arr)
filename_arr = ['hp_no_loop_'+name.split('_')[1] for name in obj2_arr]
xlabel_arr = ['Intake']*len(obj1_arr)
ylabel_arr = ['Biomass Production']*len(obj1_arr)

model_arr = ['h_pylori.json']*len(obj1_arr)


data_str = 'saturated_data/saturated_data/data'
figure_str = 'sat_figures/'

lambd = 1.0
alpha = 0.0
cutoff = 95
generations = 20
individuals = 200
nodes = 40
recon_points = 200
cores = 0
# generations = 10
# individuals = 10
# cutoff = 70

new_bounds = json.load(open('new_bounds.json'))['bounds']


for obj1_str,obj2_str,model_str,filename,xlabel,ylabel in tqdm(list(zip(obj1_arr, obj2_arr, model_arr,filename_arr,xlabel_arr,ylabel_arr))):
	###### Pareto ##########
	network = True
	model = cobra.io.load_json_model(model_str)
	for i,reaction in enumerate(model.reactions):
		reaction.bounds = new_bounds[i]
	bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))

	tqdm.write(obj2_str)

	obj1 = model.reactions.get_by_id(obj1_str).flux_expression
	if "neg" in filename:
		obj2 = -model.reactions.get_by_id(obj2_str).flux_expression
		# obj1 = -model.reactions.get_by_id(obj1_str).flux_expression
		tqdm.write('Negative')
	else:
		obj2 = model.reactions.get_by_id(obj2_str).flux_expression
		tqdm.write('Positive')
	# pops, vals, pareto_this = pareto.pareto(generations, individuals,model,obj1,obj2,batch = tqdm.write,cores = cores)

	# to_save = {'obj1_str': str(obj1), 'obj2_str': str(obj2), 'model' : model_str,'pareto': [{'obj1': p.fitness.values[0], 'obj2':p.fitness.values[1],'gene_set': list(p)} for p in pareto_this]}
	# pareto.plot_pareto(to_save['pareto'], 'pareto_'+filename, "H.Pylori Pareto Front", xlabel, ylabel, figure_str)
	# with open(data_str+filename+'.json', 'w') as outfile:
	#     json.dump(to_save, outfile)
	# tqdm.write("Pareto Generated.")

	##### Network regression #############
	try:
		network_regression.add_linear_regression(data_str+filename+'.json', cutoff)
			
	except:
		tqdm.write('Network Regression Failed')	
		network=False
	##### Pareto Reconstruction ########
	if network==True:
		# for i,x in enumerate(bounds):
		# 	model.reactions[i].bounds = x

		# pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x = pareto_reconstruction.network_reconstruct(data_str+filename+'.json', nodes, recon_points, model, obj1, obj2,cores=cores)
		# plt.clf()
		# plt.plot(pareto_x,pareto_y,'*',color='k',label='Pareto Optimal')
		# plt.plot(pareto_left[:,1], pareto_left[:,0],'r.',label='Left')
		# plt.plot(pareto_right[:,1], pareto_right[:,0],'c.',label='Right')
		# plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.',label='Noise')
		# plt.xlabel(xlabel)
		# plt.ylabel(ylabel)
		# plt.legend()
		# plt.savefig('network_recon_'+figure_str+filename+'.png')
		# plt.savefig('network_recon_'+figure_str+filename+'.eps')
		# plt.clf()

		# for i,x in enumerate(bounds):
		# 	model.reactions[i].bounds = x


		# pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x = pareto_reconstruction.reconstruct(data_str+filename+'.json', nodes, recon_points, model, obj1, obj2,cores=cores)
		# plt.clf()
		# plt.plot(pareto_x,pareto_y,'*',color='k',label='Pareto Optimal')
		# plt.plot(pareto_left[:,1], pareto_left[:,0],'r.',label='Left')
		# plt.plot(pareto_right[:,1], pareto_right[:,0],'c.',label='Right')
		# plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.',label='Noise')
		# plt.xlabel(xlabel)
		# plt.ylabel(ylabel)
		# plt.legend()
		# plt.savefig('recon_'+figure_str+filename+'.png')
		# plt.savefig('recon_'+figure_str+filename+'.eps')
		# plt.clf()

		for i,x in enumerate(bounds):
			model.reactions[i].bounds = x


		pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x = pareto_reconstruction.lin_reconstruct(data_str+filename+'.json', nodes, recon_points, model, obj1, obj2,cores=cores)
		plt.clf()
		plt.plot(pareto_x,pareto_y,'*',color='k',label='Pareto Optimal')
		plt.plot(pareto_left[:,1], pareto_left[:,0],'r.',label='Left')
		plt.plot(pareto_right[:,1], pareto_right[:,0],'c.',label='Right')
		plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.',label='Noise')
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.legend()
		plt.savefig('lin_recon_'+figure_str+filename+'.png')
		plt.savefig('lin_recon_'+figure_str+filename+'.eps')
		plt.clf()

