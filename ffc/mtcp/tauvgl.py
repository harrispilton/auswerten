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

extauin=open('m-tcp_zk.dat','r')
zk=extauin.readlines()
schmidtketau=[]
schmidtketemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	schmidtketemp.append(float(liste[0]))
	schmidtketau.append(float(liste[1]))
mytau=[]
mytemp=[]
mytauin=open('tau.dat','r')
tau=mytauin.readlines()
for i in range(0,tau.__len__()):
	liste=tau[i].split()
	mytemp.append(float(liste[0]))
	mytau.append(float(liste[1])-1)
ax.plot(mytemp,mytau,label='FFC',marker='o',linestyle='None')
ax.plot(schmidtketemp,schmidtketau,label='Schmidtke',marker='^',linestyle='None')
plt.legend()
laber=raw_input('ende')
