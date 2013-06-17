import numpy as np
from scipy.optimize import curve_fit
def func(x, a, b, c):
	return (a-x**2 + b*x + c)
	# Data (including uncertainties)
x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
y = np.array([-1.78, 4.09, 8.85, 17.9, 26.1, 35.2])
yerr = np.array([0.460, 0.658, 0.528, 1.34, 1.09, 0.786])
p0= [1.0, 3.0,-2.0]
print type(x), type(y), type(p0),type(p0[0])
# initial guess at parameters
popt, pcov= curve_fit(func, x, y, p0, yerr)
print 'optimal parameters: ', popt
print 'uncertainties of parameters: ', pcov

