import cobra
import json
import pareto_reconstruction
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm

recon_points = 200
model_name = 'h_pylori.json'

model = cobra.io.load_json_model(model_name)
# model.reactions.get_by_id('EX_o2_e').bounds = [-12.0,0.0]
	
bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))


obj1_str = 'BIOMASS_HP_published'
obj2_str = 'SUCCt2r'


obj1 = model.reactions.get_by_id(obj1_str).flux_expression
obj2 = model.reactions.get_by_id(obj2_str).flux_expression

max_bounds = np.array([-np.inf]*len(model.reactions))
min_bounds = np.array([np.inf]*len(model.reactions))
print(dir(model.optimize()))

for i,reaction in tqdm(enumerate(model.reactions),total= len(model.reactions)):
	model.objective = model.problem.Objective(reaction.flux_expression,direction='max')
	nominal_max = model.optimize()
	loopless_max = cobra.flux_analysis.loopless_solution(model)
	if loopless_max.status == 'optimal':
		min_bounds = np.minimum(loopless_max.fluxes,min_bounds)
		max_bounds = np.maximum(loopless_max.fluxes,max_bounds)
	model.objective = model.problem.Objective(reaction.flux_expression,direction='min')
	nominal_min = model.optimize()
	loopless_min = cobra.flux_analysis.loopless_solution(model)
	if loopless_min.status == 'optimal':
		min_bounds = np.minimum(loopless_min.fluxes,min_bounds)
		max_bounds = np.maximum(loopless_min.fluxes,max_bounds)
new_bounds = list(zip(min_bounds,max_bounds))

with open('new_bounds.json','w') as outfile:
	json.dump({'bounds':new_bounds}, outfile)

plt.plot(min_bounds,'.')
plt.plot(max_bounds,'.')
plt.show()
# pareto = pareto_reconstruction.origin_reconstruct(load_file, model, obj1, obj2)
# plt.plot(pareto[:,1],pareto[:,0],'r.')
# plt.xlabel(xlabel)
# plt.ylabel('Biomass Production')
# plt.title('H. Pylori Pareto Front Networked Reconstruction')
# plt.legend()
# plt.show()

# pareto_left,pareto_right,pareto_noise,pareto_y,pareto_x = pareto_reconstruction.network_reconstruct(load_file, nodes, recon_points, model, obj1, obj2,cores=cores)
# plt.figure(1)
# plt.plot(pareto_x,pareto_y,'*',color='k',label='Pareto Optimal')
# plt.plot(pareto_noise[:,1],pareto_noise[:,0],'g.',alpha = 0.8, label='Noise')
# plt.plot(pareto_right[:,1], pareto_right[:,0],'c.',alpha = 0.8, label='Right')
# plt.plot(pareto_left[:,1], pareto_left[:,0],'r.',alpha = 0.8, label='Left')
# plt.xlabel(xlabel)
# plt.ylabel('Biomass Production')
# plt.title('H. Pylori Pareto Front Networked Reconstruction')
# plt.legend()
# plt.savefig('network_recon_'+figure_save+'.png')
# plt.savefig('network_recon_'+figure_save+'.eps')

# for i,x in enumerate(bounds):
# 	model.reactions[i].bounds = x


# pareto_left_2,pareto_right_2,pareto_noise_2,pareto_y_2,pareto_x_2 = pareto_reconstruction.reconstruct(load_file, nodes, recon_points, model, obj1, obj2,cores=cores)
# plt.figure(2)
# plt.plot(pareto_x_2,pareto_y_2,'*',color='k',label='Pareto Optimal')
# plt.plot(pareto_noise_2[:,1],pareto_noise_2[:,0],'g.',alpha = 0.8, label='Noise')
# plt.plot(pareto_right_2[:,1], pareto_right_2[:,0],'c.',alpha = 0.8, label='Right')
# plt.plot(pareto_left_2[:,1], pareto_left_2[:,0],'r.',alpha = 0.8, label='Left')
# plt.xlabel(xlabel)
# plt.ylabel('Biomass Production')
# plt.title('H. Pylori Pareto Front Reconstruction')
# plt.legend()
# plt.savefig('recon_'+figure_save+'.png')
# plt.savefig('recon_'+figure_save+'.eps')
# plt.show()	



# filename_arr = ['ecoli_succ','ecoli_acet','ecoli_o2','hp_succ','hp_acet','hp_o2']
# xlabel_arr = ['Succinate Production','Acetate Production','Oxygen Intake']*2
# data_str = 'new_data/data_'
# figure_str = 'figures/pareto_'
# title_arr = ['E. Coli Pareto Front']*3
# title_arr += ['H. Pylori Pareto Front']*3

# for i,filename in enumerate(filename_arr):
# 	data = json.load(open(data_str+filename+'.json'))
# 	print(sorted(data.keys()))

# 	pareto = data['pareto']
# 	y_plot = np.array(list(map(lambda i : i['obj1'],pareto)))
# 	x_plot = np.array(list(map(lambda i : i['obj2'],pareto)))

# 	x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot,y_plot)
# 	fig = plt.figure(i)
# 	plt.xlabel(xlabel_arr[i])
# 	plt.ylabel('Biomass Production')
# 	plt.title(title_arr[i])
# 	plt.plot(x_plot,y_plot,c='b', label = 'Pareto Points',marker='.',linestyle='None')
# 	plt.plot(x0,y0,'*',c='r', label = 'Phase Transition Point')
# 	plt.plot(x_plot,kink_finder.piecewise_linear(x_plot,x0,y0,k1,k2),c='lawngreen', label = 'Fitted Line')
# 	plt.legend()
# 	plt.savefig(figure_str+filename+'.png')
# 	plt.savefig(figure_str+filename+'.eps')
# # plt.show()
