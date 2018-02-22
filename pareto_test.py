import pareto, cobra,cobra.test
import matplotlib.pyplot as plt 
import numpy as np
import json

h_p = cobra.io.load_json_model('iJO1366.json')

# print(h_p.metabolites)
# h_p = cobra.test.create_test_model("ecoli")
obj1_str = 'BIOMASS_Ec_iJO1366_WT_53p95M'
obj2_str = 'SUCCt2_2pp'
obj3_str = 'ACt4pp'
obj4_str = 'EX_o2_e'


obj1 = h_p.reactions.get_by_id(obj1_str).flux_expression
obj2 = h_p.reactions.get_by_id(obj2_str).flux_expression
obj3 = h_p.reactions.get_by_id(obj3_str).flux_expression
obj4 = h_p.reactions.get_by_id(obj4_str).flux_expression



print(h_p.reactions.get_by_id(obj1_str).gene_reaction_rule)
print(h_p.reactions.get_by_id(obj3_str).gene_reaction_rule)



print(obj1)
print(obj3)

# pareto(generations,pop_size,model, obj1, obj2, cores=0)

pops, vals, pareto = pareto.pareto(10,200,h_p,obj1,obj3)



to_save = [{'obj1': p.fitness.values[0], 'obj2':p.fitness.values[1],'gene_set': list(p)} for p in pareto]

with open('data_big_acetate.json', 'w') as outfile:
    json.dump(to_save, outfile)
print('save length')
print(len(to_save))




points = list(zip(*[p.fitness.values for p in pareto]))

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
fig.savefig('all_pareto_acetate.png')
fig.savefig('all_pareto_acetate.eps')

plt.ylabel('$f_1$')
plt.xlabel('$f_2$')
plt.title('Pareto')
plt.show()


