import cobra
import numpy as np
from matplotlib import pyplot as plt


c1 = np.array([247, 228, 251])
c2 = np.array([80, 30, 88])


def c_map(c1,c2,v_min,v_max):
	diff1 = c1-c2
	diff2 = v_max-v_min
	ratio = diff1/diff2
	return (lambda i: tuple(ratio*(i-v_min)+c2))


c_map_i = c_map(c1,c2,0,5)

print(c_map_i(1))

print(c_map_i(3))
print(c_map_i(2))

# h_p = cobra.io.load_json_model('iJO1366.json')

# print(len(h_p.metabolites))
# print(len((h_p.reactions)))
# print((len(h_p.genes)))

# print(h_p.reactions[:1000])

# print(h_p.reactions.get_by_id('BIOMASS_Ec_iJO1366_WT_53p95M').gene_reaction_rule)
# print(h_p.reactions[695])
# print(h_p.reactions[696])
# print(h_p.reactions[697])

# print(h_p.summary())
# print(h_p.genes)

# objective1 = h_p.reactions.FBA.flux_expression

# h_p.objective = h_p.problem.Objective(objective1,direction='max')
# soln = h_p.optimize()

# second = h_p.problem.Constraint(objective1,lb =soln.objective_value,ub = soln.objective_value)
# h_p.objective = h_p.problem.Objective(h_p.reactions.NH4t.flux_expression,direction='max')
# h_p.add_cons_vars(second)
# sol2 = h_p.optimize()

# print(soln.fluxes['FBA'], sol2.fluxes['FBA'])
# print(soln.fluxes['NH4t'],sol2.fluxes['NH4t'])
# BIOMASS_HP_published

# NH4t

# FBA

# FBA2

# FBP

# h_p.


# h_p = cobra.io.load_json_model('h_pylori.json')

# obj1_str = 'BIOMASS_HP_published'
# obj2_str = 'FBA'


# obj1 = h_p.reactions.get_by_id(obj1_str).flux_expression
# obj2 = h_p.reactions.get_by_id(obj2_str).flux_expression

# obj1_p = h_p.problem.Objective(obj1,direction='max')
# obj2_p = h_p.problem.Objective(obj2,direction = 'max')


# this_model = h_p.copy()
# this_model.objective = obj1_p
# # this_model.objective = this_model.problem.Objective(obj1,direction='max')
# print('obj1')

# flux1 = this_model.slim_optimize()
# print('opt1',flux1)

# this_model.objective = obj2_p
# print('obj2')

# # flux1_constr = this_model.problem.Constraint(obj1,lb=flux1,ub=flux1)
# flux1_constr = h_p.problem.Constraint(obj1,lb=flux1,ub=flux1)
# # flux1_constr = this_model.problem.Constraint(obj1,lb=flux1.objective_value,ub=flux1.objective_value)
# print('constr')


# this_model.add_cons_vars(flux1_constr)
# print('add constr')


# flux2 = this_model.slim_optimize()
# print('opt2')

