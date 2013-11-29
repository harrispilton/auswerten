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

def J_cd(omega,tau_cd,beta=1.):
	a=omega*tau_cd
	return (1/omega*(np.sin(beta*np.arctan(a)))/((1.+a**2.)**(beta/2.)))
def R1_CSA(omega,tau,delta_sigma,beta=1.):
	return 0.3*omega**2.*delta_sigma**2.*J_cd(omega,tau,beta)
def R1_DD(omega,tau,k_dd,beta=1.):
	return k_dd*(J_cd(omega,tau,beta)+4*J_cd(2*omega,tau,beta))
def R1_ges(omega,tau,k_dd,delta_sigma,beta=1.):
	return R1_DD(omega,tau,k_dd,beta)+R1_CSA(omega,tau,delta_sigma,beta)


plt.ion()

plt.figure(1)
rateax=plt.axes([0.1,0.1,0.8,0.8])
plt.xscale('log')
plt.yscale('log')
#plt.title('Einfluss von CSA bei verschiedenen Temperaturen')
rateax.set_xlabel(r'$\tau  [s]$')
rateax.set_ylabel(r'$T_1$')
plt.ylim([1e-2,1e2])
plt.xlim([2e-12,1e-7])

tau=np.logspace(-11.5,-7,120)
delta_sigma_ppm=220*1e-6

nu=121.49e6
rateax.plot(tau,1./R1_ges(nu*2*np.pi,tau/0.5,1.2e9,delta_sigma_ppm,0.5),label=r'$\nu = 121.49 MHz$',lw=2.0,color=c[0])
rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau/0.5,delta_sigma_ppm,0.5),label='CSA',color=c[0],ls='-')
rateax.plot(tau,1./R1_DD(nu*2*np.pi,tau/0.5,1.2e9,0.5),label='DD',color=c[0],ls='-.')
#rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau,delta_sigma_ppm),label='beta=1',color=c[3])

with open('jplot2/nu121.49ges.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_ges(nu*2*np.pi,t/0.5,1.2e9,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu121.49csa.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_CSA(nu*2*np.pi,t/0.5,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu121.49dd.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_DD(nu*2*np.pi,t/0.5,1.2e9,0.5))+'\n')
nu=121.49e5
rateax.plot(tau,1./R1_ges(nu*2*np.pi,tau/0.5,1.2e9,delta_sigma_ppm,0.5),label=r'$\nu = 12.149 MHz$',lw=2.0,color=c[2])
rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau/0.5,delta_sigma_ppm,0.5),label='CSA',color=c[2],ls='-')
rateax.plot(tau,1./R1_DD(nu*2*np.pi,tau/0.5,1.2e9,0.5),label='DD',color=c[2],ls='-.')
#rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau,delta_sigma_ppm),label='beta=1',color=c[3])

with open('jplot2/nu12.149ges.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_ges(nu*2*np.pi,t/0.5,1.2e9,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu12.149csa.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_CSA(nu*2*np.pi,t/0.5,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu12.149dd.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_DD(nu*2*np.pi,t/0.5,1.2e9,0.5))+'\n')

nu=121.49e7
rateax.plot(tau,1./R1_ges(nu*2*np.pi,tau/0.5,1.2e9,delta_sigma_ppm,0.5),label=r'$\nu = 1.2149 GHz$',lw=2.0,color=c[8])
rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau/0.5,delta_sigma_ppm,0.5),label='CSA',color=c[8],ls='-')
rateax.plot(tau,1./R1_DD(nu*2*np.pi,tau/0.5,1.2e9,0.5),label='DD',color=c[8],ls='-.')
#rateax.plot(tau,1./R1_CSA(nu*2*np.pi,tau,delta_sigma_ppm),label='beta=1',color=c[3])

with open('jplot2/nu1214.9ges.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_ges(nu*2*np.pi,t/0.5,1.2e9,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu1214.9csa.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_CSA(nu*2*np.pi,t/0.5,delta_sigma_ppm,0.5))+'\n')
with open('jplot2/nu1214.9dd.dat','w') as fout:
	for t in tau:
		fout.write(str(t)+' '+str(R1_DD(nu*2*np.pi,t/0.5,1.2e9,0.5))+'\n')

plt.legend()

i=raw_input('ente')
