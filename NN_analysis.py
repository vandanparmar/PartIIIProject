import numpy as np
from keras.models import load_model
from keras import backend as K
import json
from matplotlib import pyplot as plt
import pareto, cobra,cobra.test


model = load_model('NN_succinate_2.h5')

h_p = cobra.io.load_json_model('iJO1366.json')


model.summary()

data = json.load(open('data_big_succinate.json'))
X = np.array(list(map(lambda i : i['gene_set'],data)))


get_hidden_layer_output = K.function([model.layers[0].input,K.learning_phase()],[model.layers[1].output])
get_output_layer_output = K.function([model.layers[0].input,K.learning_phase()],[model.layers[3].output])
print(K.eval(model.layers[3].weights[1]))

hidden_layer_output = get_hidden_layer_output([X,0])[0]
output_layer_output = get_output_layer_output([X,0])[0]
print(hidden_layer_output.shape)

weights = K.eval(model.layers[1].weights[0]).T
print(np.shape(weights))

max_set = set()

for i,weight_set in enumerate(weights):
	if(i==0 or i==9):
		maxes = list(np.argpartition(weight_set,-10)[-10:])
		max_set.update(maxes)
		print(i,maxes, weight_set[maxes])
		print([(h_p.genes[i],h_p.genes[i].reactions) for i in maxes])
		for i in maxes:
			for reaction in h_p.genes[i].reactions:
				print(reaction)

fig, (ax1,ax2) = plt.subplots(nrows=2)

ax1.imshow(hidden_layer_output.T,aspect='auto')
ax2.imshow(output_layer_output.T,aspect='auto')
fig.savefig('NN_outputs_succ.png')
fig.savefig('NN_outputs_succ.eps')
plt.show()

print(sorted(max_set))

fig,ax = plt.subplots(nrows=1)
ax.imshow(X.T[list(max_set)],aspect='auto')
fig.savefig('express_part_succ.png')
fig.savefig('express_part_succ.eps')
plt.show()