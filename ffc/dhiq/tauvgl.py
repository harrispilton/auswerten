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
from lmfit import minimize, Parameters, report_errors
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.8,0.8])

params=Parameters()
params.add('tau0',value=-12,min=-15,max=2)
params.add('a',value=-7.0,min=-10,max=3)
params.add('einf',value=1850,min=500,max=5000)
def residuals(params,xdata,ydata=None):
	tau0=params['tau0'].value
	a=params['a'].value
	einf=params['einf'].value

	if ydata==None:
		return (tau0+((1.+np.exp(a*(-1.+(np.pi**2.*xdata)/einf)))*einf)/(xdata*np.log(10)))
	return (ydata-(tau0+((1.+np.exp(a*(-1.+(np.pi**2.*xdata)/einf)))*einf)/(xdata*np.log(10))))
	#return (ydata-(tau0+((1.+np.exp(a*(-1.+(np.pi**2.*xdata)/einf)))*einf)/(xdata*np.log(10))))

extauin=open('externe_daten/ls_TAUv1.dat','r')
zk=extauin.readlines()
ls1tau=[]
ls1temp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	ls1temp.append(float(liste[0]))
	ls1tau.append(float(liste[1]))
extauin=open('externe_daten/ls_TAUv2.dat','r')
zk=extauin.readlines()
ls2tau=[]
ls2temp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	ls2temp.append(float(liste[0]))
	ls2tau.append(float(liste[1]))
extauin=open('externe_daten/ds_tau.dat','r')
zk=extauin.readlines()
dstau=[]
dstemp=[]
for i in range(0,zk.__len__()):
	liste=zk[i].split()
	dstemp.append(float(liste[0]))
	dstau.append(float(liste[1]))


mytau=[]
mytemp=[]
mytauin=open('tau.dat','r')
tau=mytauin.readlines()
for i in range(0,tau.__len__()):
	liste=tau[i].split()
	mytemp.append(float(liste[0]))
	mytau.append(float(liste[1]))
ax.plot(mytemp,mytau,label='FFC',marker='o',linestyle='None')
###die eigenen taus sind gemalt, jetzt fuehre alle taus und temps zusammen und fitte diese liste
#for i in range(0,dstau.__len__()):
#	mytau.append(dstau[i])
#	mytemp.append(dstemp[i])
for i in range(0,ls1temp.__len__()):
	mytau.append(ls1tau[i])
	mytemp.append(ls1temp[i])
for i in range(0,ls2temp.__len__()):
	mytau.append(ls2tau[i])
	mytemp.append(ls2temp[i])
minimize(residuals,params,args=(np.array(mytemp),np.array(mytau)))
fit=residuals(params,np.array(sorted(mytemp)))
print report_errors(params)
ax.plot(sorted(mytemp),fit,label='fit')
ax.plot(ls1temp,ls1tau,label='LS1',marker='^',linestyle='None')
ax.plot(ls2temp,ls2tau,label='LS2',marker='v',ls='None')
ax.plot(dstemp,dstau,label='DS',marker='x',linestyle='None')
temp=np.linspace(200,400,100)
tau=-13.24+(980.*(1.+np.exp(-6.2*(-1.+0.004269*temp))))/temp
#ax.plot(temp,tau,label='fit')
#ax.plot(245,np.log10(1.9e-7*0.61),marker='o',ms=8,label='testpunkt')
plt.legend()
plt.show()
laber=raw_input('ende')
