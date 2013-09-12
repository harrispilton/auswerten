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

def J_debye(nu,tau):
	omega=nu*2.*np.pi
	return (tau/(1+(omega*tau)**2))
def J_cd(nu,tau,beta=1):
	omega=nu*2.*np.pi
	a=omega*tau*beta
	return (1/omega*
		np.sin(beta*np.arctan(a))/((1.+a**2)**(beta/2.))
		)
def R1_CSA(nu,tau,delta_sigma,beta=1):
	omega=nu*2.*np.pi
	return 0.3*omega**2*delta_sigma**2*J_cd(nu,tau,1)
def R1_DD(nu,tau,K_dd,beta=1):
	return K_dd*(J_cd(nu,tau,beta)+4*J_cd(2*nu,tau,beta))
def R1_ges(nu,tau,K_dd,delta_sigma,beta=1):
	return R1_DD(nu,tau,K_dd,beta)+R1_CSA(nu,tau,delta_sigma,beta)
delta_sigma_ppm=220e-6
nu=np.logspace(1.5,12.2,120)
nu=121.49e6
plt.ion()
plt.figure(1)
rateax=plt.axes([0.05,0.1,0.8,0.8])
plt.xscale('log')
plt.yscale('log')
tau=np.logspace(-11,-7,120)
rateax.plot(tau,1/R1_CSA(nu,tau,delta_sigma_ppm,0.4),label='nu ='+str(nu),color=c[1])
#rateax.plot(tau,1/R1_CSA(nu,tau,delta_sigma_ppm,0.9),color=c[2])
rateax.legend()
plt.show()
i=raw_input('ente')
