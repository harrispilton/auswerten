#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.optimize import brentq 
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
import matplotlib.cm as cm
from scipy.optimize import leastsq
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.8,0.8])

plt.yscale('log')

with open('D.dat','r') as fin:
	lines=fin.readlines()
	ds=[]
	dtemp=[]
	for line in lines:
		liste=line.split()
		dtemp.append(float(liste[0]))
		ds.append(float(liste[1]))
with open('tau.dat','r') as fin:
	lines=fin.readlines()
	tau=[]
	tautemp=[]
	for line in lines: 
		liste=line.split()
		tautemp.append(float(liste[0]))
		tau.append(float(liste[1]))
dmaltau=[]
temp=[]
for (t1,t2) in zip(dtemp,tautemp):
	print str(t1)+'    '+str(t2)
for (ta, taut, d) in zip(tau,tautemp,ds):
		dmaltau.append(10**ta*10**d/taut)
	temp.append(taut)
ax.plot(tautemp,dmaltau,marker='o',ls='None')
i=raw_input('ente')
