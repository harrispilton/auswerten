#!/usr/bin/python
#import glob
#import re
import itertools
import numpy as np
#import scipy as sp
#from scipy.optimize import curve_fit
import matplotlib.cm as cm
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider, Button, CheckButtons
c=[]
for i in np.arange(10):c.append(cm.jet(1-i/10.))
colors =itertools.cycle(c)
def Lorentz(x,f0,fwhm,a0):
	return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))

def J_debye(omega,tau):
	return (tau/(1+(omega*tau)))
def J_cd(omega,tau_cd,beta=1):
	a=omega*tau_cd
	return (np.sin(beta*np.arctan(a))/(omega*(1+a**2)**(beta/2.)))
def R1_CSA(omega,tau,delta_sigma,beta=1):
	return 0.3*omega**2*delta_sigma**2*J_cd(omega,beta)
def R1_DD(omega,tau,K_dd,beta=1):
	return K_dd*(J_cd(omega,tau,beta)+4*J_cd(2*omega,tau,beta))
def R1_ges(omega,tau,K_dd,delta_sigma,beta=1):
	return R1_DD(omega,tau,K_dd,beta)+R1_CSA(omega,tau,delta_sigma,beta)
def Susz(omega):
	return (omega*J_cd(omega,1,0.5)*(omega+40)**0.28\
		+4.e+4*(omega)*J_cd(omega+3e5,1.,0.8))\
		+1.e+9*omega*J_cd(omega+1e9,2.,1.)
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.8,0.8])
plt.xscale('log')
plt.yscale('log')
omega=np.logspace(-3,11,400)
omega=np.logspace(-3,3,300)
ax.plot(omega,omega**0.2+omega**-0.3)
ax.plot(omega,omega**0.2+2*omega**-0.3)
ax.plot(omega,omega**0.2)
#ax.plot(omega,1e4*omega**-1.2)

plt.draw()
i=raw_input('ente')
