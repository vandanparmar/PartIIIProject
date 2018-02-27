from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np
import json

def piecewise_linear(x, x0, y0, k1, k2):
	return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def get_kink_point(x,y):
	grad = (y[0]-y[-1])/(x[0]-x[-1])
	p,e = optimize.curve_fit(f=piecewise_linear,xdata = x,ydata = y,p0=[np.mean(x),np.mean(y),grad,grad])
	return p #x0,y0,k1,k2

if __name__ == '__main__':
	data = json.load(open('data_hp_succ.json'))

	y_plot = np.array(list(map(lambda i : i['obj1'],data)))
	x_plot = np.array(list(map(lambda i : i['obj2'],data)))


	x0,y0,k1,k2 = get_kink_point(x_plot,y_plot)
	print(x_plot.dtype)
	plt.plot(x_plot,y_plot,'.')
	plt.plot(x0,y0,'*',c='k')
	plt.plot(x_plot,piecewise_linear(x_plot,x0,y0,k1,k2))
	plt.show()
