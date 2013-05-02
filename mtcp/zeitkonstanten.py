#!/usr/bin/python
import glob
import re
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

plt.figure(1)
tauax=plt.axes([0.1,0.2,0.7,0.7])

axcolor = 'lightgoldenrodyellow'
axtau = plt.axes([0.8,0.1,0.4,0.02],axisbg=axcolor)
stau=Slider(axtau,'y skalierung',-2,5,valinit=1)
fzk=open('m-tcp_zk.dat','r')
dataschmidtke=fzk.readlines()
T_s=[]
tau_rot=[]
for data in dataschmidtke:
	liste = data.split()
	T_s.append(float(liste[1]))
	tau_rot.append(float(liste[0]))
plt.axes([0.1,0.2,0.7,0.7])
plt.plot(tau_rot,T_s,ls='None',marker='o')
plt.show()
