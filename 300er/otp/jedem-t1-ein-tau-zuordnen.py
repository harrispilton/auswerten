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
def J_debye(omega,tau):
	return (tau/(1+(omega*tau)))
def J_cd(omega,tau_cd,beta=1):
	a=omega*tau_cd
	return (np.sin(beta*np.arctan(a))/(omega*(1+a**2)**(beta/2.)))
def R1_DD(omega,tau,K_dd,beta=1.):
	return K_dd*(J_cd(omega,tau,beta)+4*J_cd(2*omega,tau,beta))
def R1_Q(omega,tau,K_q,beta=1.):
	return 3.*np.pi**2./10.*K_q**2.*(J_cd(omega,tau,beta)+4.*J_cd(2*omega,tau,beta))
def R2_Q(omega,tau,K_q,beta=1.):
	return 1.5*np.pi**2./10.*K_q**2.*(3.*tau+5.*J_cd(omega,tau,beta)+2.*J_cd(2.*omega,tau,beta))
def R1_ges(omega,tau,K_dd,delta_sigma,beta=1):
	return R1_DD(omega,tau,K_dd,beta)+R1_CSA(omega,tau,delta_sigma,beta)
def Susz(omega,tau,K_dd,delta_sigma,beta=1):
	return K_dd*(omega*J_cd(omega,tau,beta)+4*omega*J_cd(2*omega,tau,beta))+0.3*omega**3*delta_sigma**2*J_cd(omega,tau,beta)
delta_sigma_ppm=220e-6# dies ist eine konstante groesse angegeben in ppm diese zahl bezieht sich auf nu = 47.262 MHz und kommt aus der Diplomarbeit von Daniel Bock Seite 28
nu=np.logspace(1.5,12.2,120)
nu=[121.49e6,40e6]
K_DDs=np.logspace(1,10,5)
plt.ion()
plt.figure(1)
rateax=plt.axes([0.15,0.1,0.80,0.8])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\tau [s]$')
plt.ylabel(r'$T_1$')
#alpha=np.linspace(0,1,10)
tau=taurot=np.logspace(-10,-0,512)
#x=np.linspace(0.01,10,1e3)
#ax.plot(x,susz(x))
#plt.draw()
#plt.show()
#i = raw_input('next')
#for t in tau:
#	ax.plot(nu,R1_ges(nu,t,K_DDs[1],delta_sigma_ppm),label='tau ='+str(t)+'K_DD = 1e'+str(np.log10(K_DDs[1])),alpha=(np.log10(K_DDs[1])-2)/10,color=colors.next())
#t1dd=[1/R1_DD(46.072e6*2.*np.pi,t/0.5,3*np.pi**2/10*(182e3)**2,0.5) for t in tau]
#rateax.plot(tau,t1dd,label='nu = 46.072 MHz')
t1q=[1/R1_Q(46.072e6*2.*np.pi,t/0.5,182e3,0.5) for t in tau]
rateax.plot(tau,t1q,label='T1')
t2q=[1/R2_Q(46.072e6*2.*np.pi,t/0.5,182e3,0.5) for t in tau]
rateax.plot(tau,t2q,label='T2')
with open('t1-t2-von-tau.dat','w') as fout:
	for (t,t1,t2) in zip(tau,t1q,t2q):
		fout.write(str(t)+' '+str(t1)+' '+str(t2)+'\n')
#	rateax.plot(nu,R1_CSA(nu,tau[i],delta_sigma_ppm),label='csa',alpha=0.4)
#	rateax.plot(nu,R1_DD(nu,tau[i],K_DDs[1]),label='dd',alpha=0.4)
#x=np.logspace(-3.0,1.5,200,10)
#ax.plot(x,J_debye(x,1),label=r'$\tau =1$, $\beta=1$')
#ax.plot(x,J_debye(x,2),label=r'$\tau =2$, $\beta=1$')
#ax.plot(x,J_cd(x,5,0.2),label=r'$\tau =4$, $\beta=0.3$')
rateax.legend()
plt.show()

#suszax.plot(tau,Susz(nu[1]*2*np.pi,tau,K_DDs[1],delta_sigma_ppm,1),color=c[i])
i=raw_input('ente')
