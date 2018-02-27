import cvxpy
import json
import numpy as np
import networkx as nx
from kink_finder import get_kink_point
from matplotlib import pyplot as plt


def gene_co_express(data,cutoff):
	W = np.corrcoef(data.T)
	W = W-np.eye(np.shape(W)[0])
	W = np.power(W,2)
	clip_val = np.percentile(W.flatten(),cutoff)
	W = np.clip(W-clip_val,0,np.inf)
	to_add = np.ceil(W)
	W += clip_val*to_add
	return W


def S_from_W(W):
	W = W + np.eye(np.shape(W)[0])
	D = np.sum(W,axis=1)
	D = np.power(D,-0.5)
	D = np.diag(D)
	S = np.matmul(D,np.matmul(W,D))
	return S

data = json.load(open('data_hp_succ.json'))
ys = np.array(list(map(lambda i : i['obj1'],data)))
xs = np.array(list(map(lambda i : i['obj2'],data)))
genes = np.array(list(map(lambda i : i['gene_set'],data)))

lambd = 1.0
alpha = 0.1

W = gene_co_express(genes,95)
A = np.ceil(W)
S = S_from_W(W)

# S -= np.diag(np.diagonal(S))

def regress(genes,lambd,alpha,xs,ys,left,S):

	cost = 0
	n_genes = np.shape(genes)[1]
	constr = []

	beta = cvxpy.Variable(n_genes)
	constr.append(cvxpy.norm(beta)<=1)



	x0,y0,k1,k2 = get_kink_point(xs,ys)

	if left:
		filtered_genes = genes[ys>y0]
	else:
		filtered_genes = genes[ys<y0]


	print(x0,y0)
	print(n_genes)
	for i,(x,y,gene_set) in enumerate(zip(xs,ys,genes)):
		cost += beta.T*gene_set

	cost -= np.shape(filtered_genes)[0]*cvxpy.log_sum_exp(filtered_genes*beta)
	cost -= lambd*alpha*cvxpy.power(cvxpy.norm(beta),2)
	cost -= lambd*(1-alpha)*cvxpy.quad_form(beta,S)


	# print(cost)
	prob = cvxpy.Problem(cvxpy.Maximize(cost),constr)
	a = prob.solve(solver=cvxpy.SCS,eps=1e-5)

	print(a)
	print(beta.value)

	return beta.value





# plt.imshow(beta.value,aspect='auto')
# plt.show()

# plt.hist(beta.value,histtype='step')
# plt.show()
print(np.shape(A))
beta1 = regress(genes,lambd,alpha,xs,ys,True,S).flatten()
beta2 = regress(genes,lambd,alpha,xs,ys,False,S).flatten()


to_save = {'beta1':beta1.tolist(),'beta2':beta2.tolist(),'A':A.tolist()}
data['network'] = to_save
with open('network_hp_succ.json', 'w') as outfile:
    json.dump(data, outfile)

print(np.shape(beta1))


plt.figure(1)
maxes1 = np.argpartition(beta1,-20)[-20:]
print(np.shape(A[maxes1][:,maxes1]))
print(maxes1)
small1 = A[maxes1][:,maxes1]
G = nx.from_numpy_matrix(small1)
mapping=dict(zip(G.nodes(),maxes1))
G = nx.relabel_nodes(G,mapping)
nx.draw_spring(G,with_labels=True,font_color='k')



plt.figure(2)


maxes2 = np.argpartition(beta2,-20)[-20:]
print(maxes2)
print(np.shape(A[maxes2][:,maxes2]))
small2 = A[maxes2][:,maxes2]
G = nx.from_numpy_matrix(small2)
mapping=dict(zip(G.nodes(),maxes2))
G = nx.relabel_nodes(G,mapping)
nx.draw_spring(G,with_labels=True,font_color='k')
plt.show() 
