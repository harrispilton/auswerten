#!/usr/bin/python
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
beta=0.9
omega=np.logspace(5.0,8.5,20,10)
tau_c=1.0
K_DD=1.0
delta_sigma_CSA=0.226
sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
def J(omega,tau_c):
	return tau_c/(1. + omega **2 * tau_c **2)**beta
def Chi(omega,tau_c):
	return omega*J(omega,tau_c)
plt.axes([0.1,0.1,0.55,0.8])
sefdata=[]
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	chi=[]
	brlx=[]
	t1=[]
	r1=[]
	percerr=[]
	abserr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		brlx.append(liste[0])
		brlx=map(float,brlx)
		t1.append(liste[1])
		t1=map(float,t1)
		r1.append(liste[2])
		r1=map(float,r1)
		percerr.append(liste[3])
		percerr=map(float,percerr)
		abserr.append(liste[4])
		abserr=map(float,abserr)
		zone.append(liste[5])
		zone=map(int,zone)
		relativefile.append(liste[6])
	for i,b in enumerate(brlx):
		brlx[i]=brlx[i]*10e6
		chi.append(r1[i]*brlx[i])
	plt.plot(brlx,chi,label=relativefile[0])
print Chi(20e7,1e-9)
plt.plot(omega, Chi(omega,1e-6),label='chi mit tau_c =1')

plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('omega')
plt.ylabel('Chi')
plt.show()
