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
from lmfit import  minimize, Parameters, report_errors

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
extauin=open('externe_daten/mTCP_PCS_Tau.dat','r')
zk=extauin.readlines()
pcstau=[]
pcstemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	pcstemp.append(1000./float(liste[0]))
	pcstau.append(float(liste[1]))

mytau=[]
mytemp=[]
mytauin=open('tau.dat','r')
tau=mytauin.readlines()
for i in range(0,tau.__len__()):
	liste=tau[i].split()
	mytemp.append(float(liste[0]))
	mytau.append(float(liste[1]))
ax.plot(mytemp,mytau,label='FFC',marker='o',linestyle='None')
ax.plot(schmidtketemp,schmidtketau,label='Schmidtke',marker='^',linestyle='None')
ax.plot(pcstemp,pcstau,label='pcs')
temp=np.linspace(200,400,100)
tau=-13.24+(980.*((1.+np.exp(-6.2*(-1.+0.004269*temp)))))/temp
ax.plot(temp,tau,label='fit')
plt.legend()
laber=raw_input('ende')
