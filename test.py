import json
import kink_finder
from matplotlib import pyplot as plt
import numpy as np

filename_arr = ['ecoli_succ','ecoli_acet','ecoli_o2','hp_succ','hp_acet','hp_o2']
xlabel_arr = ['Succinate Production','Acetate Production','Oxygen Intake']*2
data_str = 'new_data/data_'
figure_str = 'figures/pareto_'
title_arr = ['E. Coli Pareto Front']*3
title_arr += ['H. Pylori Pareto Front']*3

for i,filename in enumerate(filename_arr):
	data = json.load(open(data_str+filename+'.json'))
	print(sorted(data.keys()))

	pareto = data['pareto']
	y_plot = np.array(list(map(lambda i : i['obj1'],pareto)))
	x_plot = np.array(list(map(lambda i : i['obj2'],pareto)))

	x0,y0,k1,k2 = kink_finder.get_kink_point(x_plot,y_plot)
	fig = plt.figure(i)
	plt.xlabel(xlabel_arr[i])
	plt.ylabel('Biomass Production')
	plt.title(title_arr[i])
	plt.plot(x_plot,y_plot,c='b', label = 'Pareto Points',marker='.',linestyle='None')
	plt.plot(x0,y0,'*',c='r', label = 'Phase Transition Point')
	plt.plot(x_plot,kink_finder.piecewise_linear(x_plot,x0,y0,k1,k2),c='lawngreen', label = 'Fitted Line')
	plt.legend()
	plt.savefig(figure_str+filename+'.png')
	plt.savefig(figure_str+filename+'.eps')
# plt.show()
