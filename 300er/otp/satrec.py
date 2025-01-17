#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
#import scipy as sp
#from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from lmfit import minimize, Parameters, report_errors
#import matplotlib.cm as cm
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider, Button, CheckButtons

def residual(params,xdata,ydata=None):
	t1=params['t1'].value
	t0=params['t0'].value
	satval=params['satval'].value
	if ydata is None:
		return satrec(xdata,satval,t1,t0)
	return (ydata-satrec(xdata,satval,t1,t0))
def Lorentz(x,f0,hwhm,a0,y0=0):
	return (y0+a0/(1.+((x-f0)/hwhm)**2))
	#return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))
def satrec(t,satval,t1,t0):
	return satval*(-(np.exp(-(t+t0)/t1))+1)

params = Parameters()
params.add('logt1',value=-2.0,min=-8.0,max=5.0)
pi=np.pi
params.add('t1',expr='(10**(logt1))')
params.add('t0',value=0,min=-1,max=1)
params.add('satval', value=1, min=0.85,max=1.2,vary=False)
#params.add('f0', value=0,min=-1e4,max=1e4)
#params.add('y0',value=0.05,min=0,max=0.2)
plt.ion()
x=np.linspace(-10,10,1e3)
markers = itertools.cycle(['o','s','v'])
files=glob.glob('otp/satrec/'+'*K.dat')
files.sort()
p0=[2,1.0e4,1.0]
#a0_a=[]
#fwhm_a=[]
#f0_a=[]
#for freq in fr:
#	a0_a.append(a0)
#	fwhm_a.append(fwhm)
#	f0_a.append(f0)
times=[]
inns=[]
temps=[]
t1s=[]
for data in files:
	plt.figure(1)
	ax=plt.axes([0.1,0.1,0.8,0.8])
	plt.xscale('log')
	f=open(data,'r')
	lines=f.readlines()
	for i in range(0,2): lines.pop(0)
	time=[]
	inn=[]
	for line in lines:
		liste=line.split()
		time.append(np.float(liste[0]))
		inn.append((np.float(liste[1])**2+np.float(liste[2])**2)**0.5)
	maxinn=max(inn)
	inn=[i/maxinn for i in inn]
	times.append(time)
	inns.append(inn)
	temps.append(float(data[11:-6]))
	print temps
	ax.plot(time,inn,label=str(data[11:-6]),marker=markers.next(),linestyle='None')	
	inn=np.array(inn)
	time=np.array(time)
	out = minimize(residual, params,args=(time,inn))
	fit=residual(params,np.array(sorted(time)))
	ax.plot(sorted(time),fit,label=str(temps[-1])+' Fit')
	print report_errors(params)
	plt.legend()
	plt.autoscale()
	plt.draw()
	plt.show()
	t1s.append(params['t1'].value)
	i=raw_input('next')
plt.legend()
plt.show()
with open('satrec.dat','w') as fout:
	for i in range(0,temps.__len__()):
		fout.write(str(temps[i])+' '+str(t1s[i])+'\n')	
i=raw_input('ente')
