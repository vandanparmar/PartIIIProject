from pareto_reconstruction import *
import os
import json
from tqdm import tqdm
import numpy as np
from collections import Counter
import csv

filepath = 'saturated_data/ratios/'

files = os.listdir(filepath)

writer = csv.writer(open("ratios.tsv", 'w'),delimiter='\t')
writer.writerow(['name','left','noise','right','original'])

for file in tqdm(files):
	data = json.load(open(filepath+file))
	if 'no_loop' not in file:
		continue

	# pareto_data = data['pareto']
	# pareto = np.array(list(map(lambda x: [x['obj1'],x['obj2']],pareto_data)))


	# if 'network_recon' in data:
	# 	paretos = prep_fba_set(np.array(data['network_recon']['pareto_left']), 'left')
	# 	paretos.extend(prep_fba_set(np.array(data['network_recon']['pareto_right']),'right'))
	# 	paretos.extend(prep_fba_set(np.array(data['network_recon']['pareto_noise']),'noise'))
	# 	# paretos.extend(prep_fba_set(pareto,'original'))
	# 	pareto_collection = find_optimal(paretos)
	# 	parts = list(map(lambda x: x[2],pareto_collection))
	# 	ratios = dict(Counter(parts))
	# 	data['network_recon']['ratios'] = ratios
	# if 'recon' in data:
	# 	paretos = prep_fba_set(np.array(data['recon']['pareto_left']), 'left')
	# 	paretos.extend(prep_fba_set(np.array(data['recon']['pareto_right']),'right'))
	# 	paretos.extend(prep_fba_set(np.array(data['recon']['pareto_noise']),'noise'))
	# 	# paretos.extend(prep_fba_set(pareto,'original'))
	# 	pareto_collection = find_optimal(paretos)
	# 	parts = list(map(lambda x: x[2],pareto_collection))
	# 	ratios = dict(Counter(parts))
	# 	data['recon']['ratios'] = ratios
	# if 'lin_recon' in data:
	# 	paretos = prep_fba_set(np.array(data['lin_recon']['pareto_left']), 'left')
	# 	paretos.extend(prep_fba_set(np.array(data['lin_recon']['pareto_right']),'right'))
	# 	paretos.extend(prep_fba_set(np.array(data['lin_recon']['pareto_noise']),'noise'))
	# 	# paretos.extend(prep_fba_set(pareto,'original'))
	# 	pareto_collection = find_optimal(paretos)
	# 	parts = list(map(lambda x: x[2],pareto_collection))
	# 	ratios = dict(Counter(parts))
	# 	data['lin_recon']['ratios'] = ratios

	# with open('saturated_data/ratios/'+file,'w') as f:
	# 	json.dump(data, f)


	if 'network_recon' in data:
		d = data['network_recon']['ratios']
		w = [0]*5
		w[0] = file+' network'
		if 'left' in d:
			w[1] = d['left']
		if 'noise' in d:
			w[2] = d['noise']
		if 'right' in d:
			w[3] = d['right']
		if 'original' in d:
			w[4] = d['original']
		writer.writerow(w)

	if 'recon' in data:
		d = data['recon']['ratios']
		w = [0]*5
		w[0] = file+' recon'
		if 'left' in d:
			w[1] = d['left']
		if 'noise' in d:
			w[2] = d['noise']
		if 'right' in d:
			w[3] = d['right']
		if 'original' in d:
			w[4] = d['original']
		writer.writerow(w)

	if 'lin_recon' in data:
		d = data['lin_recon']['ratios']
		w = [0]*5
		w[0] = file+' linear'
		if 'left' in d:
			w[1] = d['left']
		if 'noise' in d:
			w[2] = d['noise']
		if 'right' in d:
			w[3] = d['right']
		if 'original' in d:
			w[4] = d['original']
		writer.writerow(w)
