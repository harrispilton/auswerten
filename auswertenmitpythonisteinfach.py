#!/usr/bin/python
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
beta=0.6
print 'hallo'
#omega=1.0
omega=np.logspace(4.0,9.0,300,10)
tau_c=1e-30
K_DD=1.0
delta_sigma_CSA=0.226
def K_CSA(omega):
	a=2.0/18*(omega * delta_sigma_CSA)**2
	return a
def J(omega,tau_c):
	a=tau_c * 1.0/(1.0+(omega * tau_c)**2)**beta
	return a
def R_1(omega,tau_c):
	a=1/(K_DD * J(omega,tau_c) + K_CSA(omega) * J(omega,tau_c))
	return a
alle=[]
temp=glob.glob('*.dat')
temp.remove('all.dat')
temp.sort()
brlxs=[]
chis=[]
for i, val in enumerate(temp):
	x, y1, y2 = np.loadtxt(temp[i], unpack=True, usecols=[0,1,2])
	brlx = list(x)
	t1 = list(y1)
	r1 = list(y2)
	cbrlx=[]
	ct1=[]
	cr1=[]
	chi=[]
	for ii, val in enumerate(brlx):
		if t1[ii]>0.002:
			cbrlx.append(brlx[ii])
			ct1.append(t1[ii])
			cr1.append(r1[ii])
			chi.append(0)
	for ii, val in enumerate(cbrlx):
		cbrlx[ii]=10e6*cbrlx[ii]
		chi[ii]=10e6*cbrlx[ii]/ct1[ii]	
	brlxs.append(cbrlx)
	chis.append(chi)
#if t1s in t1 > 0.001: plt.plot(brlx,t1)
	plt.plot(cbrlx, chi,label=temp[i])
plt.plot(omega, R_1,'r--')
plt.legend(loc='upper left')
plt.yscale('log')
plt.xscale('log')
plt.show()

for i,val in enumerate(temp):
	f = open(temp[i],'r')
	lines=f.readlines()
#	alle.append[lines],
