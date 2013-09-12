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
	tau_cd=tau/beta
	return 0.3*omega**2*delta_sigma*J_cd(omega,tau_cd,beta)
def R1_DD(omega,tau,K_dd,beta=1):
	return K_dd*(J_cd(omega,tau_cd,beta)+4*J_cd(2*omega,tau_cd,beta))
def R1_ges(omega,tau,K_dd,delta_sigma,beta=1):
	return R1_DD(omega,tau,K_dd,beta)+R1_CSA(omega,tau,delta_sigma,beta)
def Susz(omega,tau,K_dd,delta_sigma,beta=1):
	tau_cd=tau/beta
	return K_dd*(omega*J_cd(omega,tau_cd,beta)+4*omega*J_cd(2*nu,tau_cd,beta))+0.3*omega**3*delta_sigma*J_cd(omega,tau_cd,beta)
delta_sigma_ppm=220e-6# dies ist eine konstante groesse angegeben in ppm diese zahl bezieht sich auf nu = 47.262 MHz und kommt aus der Diplomarbeit von Daniel Bock Seite 28
nu=np.logspace(1.5,12.2,120)
nu=[121.49e9,40e9]
K_DDs=np.logspace(1,10,5)
plt.ion()
plt.figure(1)
rateax=plt.axes([0.05,0.1,0.35,0.8])
plt.xscale('log')
plt.yscale('log')

suszax=plt.axes([0.5,0.1,0.4,0.8])
plt.xscale('log')
plt.yscale('log')
#alpha=np.linspace(0,1,10)
tau=np.logspace(-11,-7,120)
#x=np.linspace(0.01,10,1e3)
#ax.plot(x,susz(x))
#plt.draw()
#plt.show()
#i = raw_input('next')
#for t in tau:
#	ax.plot(nu,R1_ges(nu,t,K_DDs[1],delta_sigma_ppm),label='tau ='+str(t)+'K_DD = 1e'+str(np.log10(K_DDs[1])),alpha=(np.log10(K_DDs[1])-2)/10,color=colors.next())

for i in range(0,nu.__len__()):
	rateax.plot(tau,1/R1_ges(nu[i],tau,K_DDs[1],delta_sigma_ppm,0.4),label='nu ='+str(nu[i])+'K_DD = 1e'+str(np.log10(K_DDs[4])),
	alpha=1,#(np.log10(K_DDs[4])-2)/10,
	color=c[i])
	rateax.plot(tau,1/R1_ges(nu[i],tau,K_DDs[1],delta_sigma_ppm,0.9),label='nu ='+str(nu[i])+'K_DD = 1e'+str(np.log10(K_DDs[4])),
	alpha=1,#(np.log10(K_DDs[4])-2)/10,
	color=c[i])

	suszax.plot(nu,Susz(nu,tau[i],K_DDs[1],delta_sigma_ppm,1),color=c[i])
#	rateax.plot(nu,R1_CSA(nu,tau[i],delta_sigma_ppm),label='csa',alpha=0.4)
#	rateax.plot(nu,R1_DD(nu,tau[i],K_DDs[1]),label='dd',alpha=0.4)
#x=np.logspace(-3.0,1.5,200,10)
#ax.plot(x,J_debye(x,1),label=r'$\tau =1$, $\beta=1$')
#ax.plot(x,J_debye(x,2),label=r'$\tau =2$, $\beta=1$')
#ax.plot(x,J_cd(x,5,0.2),label=r'$\tau =4$, $\beta=0.3$')
rateax.legend()
plt.show()

i=raw_input('ente')
