import cobra
import json
import pareto_reconstruction
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
from cobra.flux_analysis import flux_variability_analysis
model_name = 'h_pylori.json'

model = cobra.io.load_json_model(model_name)
bounds = flux_variability_analysis(model,model.reactions).as_matrix()

bounds = np.array(list(map(lambda reaction : reaction.bounds, model.reactions)))
# bounds[:,1] -= 10
# bounds[:,0] += 10
# bounds[:,1] = np.clip(bounds[:,1],-np.inf,np.inf)
# bounds[:,0] = np.clip(bounds[:,0],-np.inf,np.inf)

# for i,x in enumerate(bounds):
# 	cur_bound = model.reactions[i].bounds
# 	if cur_bound[0]==0.0:
# 		bounds[i,1] = 0.0
# 		x[1] = 0.0
# 	if cur_bound[1]==0.0:
# 		bounds[i,0] = 0.0
# 		x[0] = 0.0

# 	model.reactions[i].bounds = x[[1,0]]
	
# 	if bounds[i,0]<0:
# 		bound = bounds[i,0]
# 		bounds[i,0] = bounds[i,1]
# 		bounds[i,1] = bound

# print(bounds)

for i,reaction in enumerate(model.reactions):
	if 'EX_' in reaction.id:
		print(reaction.id)
		print(reaction.bounds)
		print(bounds[i][[1,0]])
		reaction.bounds = np.clip(reaction.bounds, -12.0, 12.0)

unclipped = {}
for reaction in tqdm(model.reactions):
	# if reaction.id=='ACKr':
	# 	print(reaction.bounds)
	model.objective = model.problem.Objective(reaction.flux_expression,direction='max')
	soln = model.optimize()
	# print(model.reactions.get_by_id('ACKr').flux)

	upper = soln.fluxes/bounds[:,0]
	upper[upper==np.inf] = 0
	upper[upper==1.0] = 0
	max_upper = np.max(upper)
	upper_reac = np.where(max_upper==upper)[0]
	if len(upper_reac)==len(model.reactions):
		upper_reac = []
	upper_names = list(map(lambda i : model.reactions[int(i)].id,upper_reac))


	lower = soln.fluxes/bounds[:,1]
	lower[lower==np.inf] = 0
	lower[lower==1.0] = 0
	max_lower = np.max(lower)
	lower_reac = np.where(max_lower==lower)[0]
	if len(lower_reac)==len(model.reactions):
		lower_reac = []
	lower_names = list(map(lambda i : model.reactions[int(i)].id,lower_reac))

	unclipped[reaction.id] = {'bounds':reaction.bounds,'val':reaction.flux,'upper_val':max_upper,'upper_names':upper_names,'lower_val':max_lower,'lower_names':lower_names}

print(unclipped)
with open('e_reac_saturated.json','w') as outfile:
	json.dump(unclipped, outfile) 