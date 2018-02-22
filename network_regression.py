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

data = json.load(open('data_big_succinate.json'))
ys = np.array(list(map(lambda i : i['obj1'],data)))
xs = np.array(list(map(lambda i : i['obj2'],data)))
genes = np.array(list(map(lambda i : i['gene_set'],data)))


W = gene_co_express(genes,99)




A = np.ceil(W)

S = S_from_W(W)

# S -= np.diag(np.diagonal(S))

cost = 0

n_genes = np.shape(genes)[1]

lambd = 1.0
alpha = 0.0

beta = cvxpy.Variable(n_genes)
constr = [cvxpy.norm(beta)<=1]



x0,y0,k1,k2 = get_kink_point(xs,ys)

filtered_genes = genes[ys>y0]

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

plt.imshow(beta.value,aspect='auto')
plt.show()

plt.hist(beta.value,histtype='step')
plt.show()