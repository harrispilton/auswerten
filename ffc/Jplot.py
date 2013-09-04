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
	return (tau/(1+(omega*tau)))
def J_cd(nu,tau,beta=0.5):
	omega=nu*2.*np.pi
	a=omega*beta*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1+a)**(beta/2)))
def R1_CSA(nu,tau,delta_sigma):
	omega=nu*2.*np.pi
	return 0.3*omega**2*J_cd(nu,tau)
def R1_DD(nu,tau,K_dd):
	return K_dd*(J_cd(nu,tau)+4*J_cd(2*nu,tau))
def R1_ges(nu,tau,K_dd,delta_sigma):
	return R1_DD(nu,tau,K_dd)+R1_CSA(nu,tau,delta_sigma)
delta_sigma_ppm=226# dies ist eine fast konstante groesse angegeben in ppm diese zahl bezieht sich auf nu = 47.262 MHz und kommt aus der Diplomarbeit von Daniel Bock Seite 28
#umrechnung fuer unser system:
nu=np.logspace(1.5,10.2,120)
K_DDs=np.logspace(6,10,5)
plt.figure(1)
ax=plt.axes([0.1,0.15,0.8,0.8])
plt.xscale('log')
plt.yscale('log')
alpha=np.linspace(0,1,10)
tau=np.logspace(-11,-7,5)
print np.log10(K_DDs)
for t in tau:
	ax.plot(nu,R1_ges(nu,t,K_DDs[1],delta_sigma_ppm*1e-6/(nu*2*np.pi)),label='tau ='+str(t)+'K_DD = 1e'+str(np.log10(K_DDs[1])),alpha=(np.log10(K_DDs[1])-2)/10,color=colors.next())

for i in range(0,tau.__len__()):
	ax.plot(nu,R1_ges(nu,tau[i],K_DDs[4],delta_sigma_ppm*1e-6/(nu*2*np.pi)),label='tau ='+str(tau[i])+'K_DD = 1e'+str(np.log10(K_DDs[4])),alpha=(np.log10(K_DDs[4])-2)/10,color=c[i])
	ax.plot(nu,R1_CSA(nu,tau[i],delta_sigma_ppm*1e-6/(nu*2*np.pi)),label='csa',alpha=0.4)
	ax.plot(nu,R1_DD(nu,tau[i],K_DDs[4]),label='dd',alpha=0.4)
#x=np.logspace(-3.0,1.5,200,10)
#ax.plot(x,J_debye(x,1),label=r'$\tau =1$, $\beta=1$')
#ax.plot(x,J_debye(x,2),label=r'$\tau =2$, $\beta=1$')
#ax.plot(x,J_cd(x,5,0.2),label=r'$\tau =4$, $\beta=0.3$')
plt.legend()
plt.show()


