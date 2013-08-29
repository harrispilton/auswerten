#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.optimize import brentq 
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.optimize import leastsq
from lmfit import  minimize, Parameters, report_errors

p_true = Parameters()
p_true.add('amp', value=14.0)
p_true.add('period', value=5.33)
p_true.add('shift', value=0.123)
p_true.add('decay', value=0.010)

def residual(pars, x, data=None):
	amp = pars['amp'].value
	per = pars['period'].value
	shift = pars['shift'].value
	decay = pars['decay'].value
	
	if abs(shift) > np.pi/2:
		shift = shift - np.sign(shift)*np.pi
	model = amp*np.sin(shift + x/per) * np.exp(-x*x*decay*decay)
	if data is None:
		return model
	return (model - data)

n = 20
xmin = 10.
xmax = 250.0
noise = np.random.normal(scale=0.7215, size=n)
x = np.linspace(xmin, xmax, n)
data  = residual(p_true, x) + noise

fit_params = Parameters()
fit_params.add('amp', value=13.0)
fit_params.add('period', value=2)
fit_params.add('shift', value=0.0)
fit_params.add('decay', value=0.02)
out = minimize(residual, fit_params, args=(x,), kws={'data':data})
fit = residual(fit_params, x)
report_errors(fit_params)
plt.plot(x,fit)
plt.plot(x,data,linestyle='None',marker='o')
plt.show()
