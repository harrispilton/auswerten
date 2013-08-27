#!/usr/bin/python
#import glob
#import re
#import itertools
import numpy as np
#import scipy as sp
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider, Button, CheckButtons

def Lorentz(x,f0,fwhm,a0):
	return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))

def J_debye(omega,tau):
    return (tau/(1+(omega*tau)))
def J_cd(omega,tau,beta):
    a=omega*beta*tau
    return (np.sin(beta*np.arctan(a))/(omega*(1+a)**(beta/2)))
plt.figure(1)
ax=plt.axes([0.1,0.15,0.8,0.8])
#ax.logscale()
x=np.logspace(-3.0,1.5,200,10)
ax.plot(x,J_debye(x,1),label=r'$\tau =1$, $\beta=1$')
ax.plot(x,J_debye(x,2),label=r'$\tau =2$, $\beta=1$')
ax.plot(x,J_cd(x,5,0.2),label=r'$\tau =4$, $\beta=0.3$')
plt.legend()
plt.show()


