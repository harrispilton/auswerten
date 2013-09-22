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

def residual(params,xdata,ydata=None):
	p1=params['p1'].value
	p2=params['p2'].value
	p3=params['p3'].value
	p4=params['p4'].value
	if ydata is None: 
		return (p1+(p2*(1.+np.exp(p3*(-1.+p4*xdata))))/xdata)
	return (ydata-(p1+(p2*(1.+np.exp(p3*(-1.+p4*xdata))))/xdata))
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.8,0.8])
alltau=[]
alltemp=[]
extauin=open('externe_daten/OTP - Decalin 85 zu 15_LS_ZK.txt','r')
zk=extauin.readlines()
schmidtketau=[]
schmidtketemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	schmidtketemp.append(float(liste[0]))
	schmidtketau.append(float(liste[1]))


extauin=open('externe_daten/otp_NikoDLS.dat','r')

zk=extauin.readlines()
pcstau=[]

pcstemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	pcstemp.append(float(liste[0]))
	pcstau.append(float(liste[1]))
extauin.close()
extauin=open('externe_daten/OTP_PCS.dat','r')
zk=extauin.readlines()
dlstau=[]
dlstemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	dlstau.append(float(liste[1]))
	dlstemp.append(float(liste[0]))
params=Parameters()
params.add('p1',value = -12,min=-20,max=0)
params.add('p2',value = 930,min=10,max=10000)
params.add('p3',value=-3,min=-10,max=-1)
params.add('p4',value=0.04,min=0.000001,max=1)



mytau=[]
mytemp=[]
mytauin=open('tau.dat','r')
tau=mytauin.readlines()
for i in range(0,tau.__len__()):
	liste=tau[i].split()
	mytemp.append(float(liste[0]))
	mytau.append(float(liste[1]))
ax.plot(mytemp,mytau,label='FFC, Decalin in 85%otp',marker='o',linestyle='None')
ax.plot(schmidtketemp,schmidtketau,label='LS 85%otp',marker='^',linestyle='None')
ax.plot(pcstemp,pcstau,label='PCS reines otp',marker='v',ls='None')
ax.plot(dlstemp,dlstau,label='DLS reines otp',marker='x',ls='None')
temp=np.linspace(200,400,100)
tau=-12.24+(980.*(1.+np.exp(-6.2*(-1.+0.004269*temp))))/temp
out=minimize(residual, params,args=(np.array(mytemp),np.array(mytau)))
print report_errors(params)
fit=residual(params,np.array(mytemp))
ax.plot(mytemp,fit,label='fit')
#ax.plot(245,np.log10(1.9e-7*0.61),marker='o',ms=8,label='testpunkt')
plt.legend()
laber=raw_input('ende')
