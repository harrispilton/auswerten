#!/usr/bin/python
import glob
import re
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

def yscale(faktor):
	D=[]
	for ds in d:
		D.append(ds+faktor)
	tauax.lines[1].set_ydata(D)

plt.figure(1)
tauax=plt.axes([0.1,0.2,0.7,0.7])

axcolor = 'lightgoldenrodyellow'
axtau = plt.axes([0.2,0.1,0.4,0.02],axisbg=axcolor)
stau=Slider(axtau,'y skalierung',-30,-10,valinit=1)
stau.on_changed(yscale)
fzk=open('m-tcp_zk.dat','r')
dataschmidtke=fzk.readlines()
fD=open('m-tcp_D.dat','r')
dif=fD.readlines()
T_s=[]
tau_rot=[]
for data in dataschmidtke:
	liste = data.split()
	T_s.append(float(liste[0]))
	tau_rot.append((float(liste[1])))
print T_s
T_d=[]
d=[]
for data in dif:
	liste = data.split()
	T_d.append(float(liste[0]))
	d.append(-(float(liste[1])))
plt.axes([0.1,0.2,0.7,0.7])
plt.plot(T_s,tau_rot,ls='None',marker='x',label='Schmidtke Daten')
plt.plot(T_d,d,ls='None',marker='o',label='Daten FFC (skaliert)')
plt.legend(loc='upper right')
plt.show()
