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
	f0=params['f0'].value
	hwhm=10**params['loghwhm'].value
	a0=params['a0'].value
	y0=params['y0'].value
	if ydata is None:
		return Lorentz(xdata,f0,hwhm,a0,y0)
	return (ydata-Lorentz(xdata,f0,hwhm,a0,y0))
def Lorentz(x,f0,hwhm,a0,y0=0):
	return (y0+a0/(1.+((x-f0)/hwhm)**2))
	#return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))
params = Parameters()
params.add('loghwhm',value=2.0,min=-20.0,max=20.0)
pi=np.pi
params.add('t_2',expr='(2/(pi*10**(loghwhm)))')
params.add('a0', value=1, min=0.95,max=1.1,vary=False)
params.add('f0', value=0,min=-1e4,max=1e4)
params.add('y0',value=0.05,min=0,max=0.2)
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.8,0.8])
x=np.linspace(-10,10,1e3)
markers = itertools.cycle(['o','s','v'])
files=glob.glob('otp/1d/'+'*K.dat')
files.sort()
print type(x)
p0=[2,1.0e4,1.0]
#a0_a=[]
#fwhm_a=[]
#f0_a=[]
#for freq in fr:
#	a0_a.append(a0)
#	fwhm_a.append(fwhm)
#	f0_a.append(f0)
freqs=[]
betrags=[]
temps=[]
fittemps=[]
t2s=[]
for data in files:
	f=open(data,'r')
	lines=f.readlines()
	for i in range(0,2): lines.pop(0)
	freq=[]
	betrag=[]
	for line in lines:
		liste=line.split()
		freq.append(np.float(liste[0]))
		betrag.append(((np.float(liste[1]))**2+(np.float(liste[2])**2))**0.5)
	maxbetrag=max(betrag)
	betrag=[b/maxbetrag for b in betrag]
	freqs.append(freq)
	betrags.append(betrag)
	temps.append(float(data[7:-6]))
	print temps
	ax.plot(freq,betrag,label=str(data),marker=markers.next(),linestyle='None')	
	plt.autoscale()
	plt.draw()
	#plt.plot(fr,Lorentz(fr,f0,fwhm,a0))
	fitten='y'
	if str(fitten)=='y':
		params['loghwhm'].value=2
		#for i in range(freq.__len__()-1,-1,-1):
		#	if abs(freq[i])>1.5e4:
		#		freq.pop(i)
		#		betrag.pop(i)
		freq=np.array(freq)
		betrag=np.array(betrag)
		out = minimize(residual, params,args=(freq,betrag))
		result=freq+out.residual
		fit = residual(params, np.array(freq))
		print report_errors(params)
		ax.plot(freq,fit,label=str(temps[-1])+' Fit')
		plt.legend()
		t2s.append(params['t_2'].value)
		fittemps.append(temps[-1])

with open('otp_T2.dat','w') as fout: 
	for i in range(0,t2s.__len__()):
		fout.write(str(fittemps[i])+' '+str(t2s[i])+'\n')
#	print type(freq), type(betrag), type(p0), type(p0[0])
#	fitpars, covmat = curve_fit(
#			Lorentz,
#			freq,
#			betrag,
#			p0)
#			#[f0,fwhm,a0],
#			#maxfev=10000
#			#)
#	print 'fitparameter '+str(fitpars)+'\n\ncovmat '+str(covmat)+'\n\n\n'
#	print "fitpars[0]"+str(fitpars[0])
#	print "fitpars[0]"+str(fitpars[0])
#	print "fitpars[0]"+str(fitpars[0])
#	print "fitpars[0]"+str(fitpars[0])
#	with open('otp_T1.dat','a') as fout:
#		fout.write(str(data[7:10])+'\t'+str(np.pi/fitpars[1])+'\n')
#	plt.plot(freq,Lorentz(freq,fitpars[0],fitpars[1],fitpars[2]),
			#label=data[7:10]
#			)
#	p0=fitpars
plt.legend()
plt.show()


plt.figure(2)
plt.xscale('linear')
plt.yscale('log')
plt.plot(fittemps,t2s,linestyle='None',marker='x')
plt.draw()
i=raw_input('ente')
