#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

def Lorentz(x,f0,fwhm,a0):
	if(fwhm<1e-6):
		return 1/x
	else:
		return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))

markers = itertools.cycle(['o','s','v'])
files=glob.glob('otp/1d/'+'*K.dat')
files.sort()
x=np.linspace(-1e6,1e6,1e5)
print type(x)
p0=[2e5,1.0e4,2.0]
#a0_a=[]
#fwhm_a=[]
#f0_a=[]
#for freq in fr:
#	a0_a.append(a0)
#	fwhm_a.append(fwhm)
#	f0_a.append(f0)

with open('otp_T1.dat','w') as fout: fout.close()
for data in files:
	f=open(data,'r')
	lines=f.readlines()
	for i in range(0,2): lines.pop(0)
	freq=[]
	betrag=[]
	for line in lines:
		liste=line.split()
		freq.append(np.float(liste[0]))
		betrag.append(((np.float(liste[1]))**2+(np.float(liste[2])**2)**0.5))
	plt.figure(1)
	plt.plot(freq,betrag,label=str(data),marker=markers.next(),linestyle='None')	
	#plt.plot(fr,Lorentz(fr,f0,fwhm,a0))
	freq=np.array(freq)
	betrag=np.array(betrag)
	print type(freq), type(betrag), type(p0), type(p0[0])
	fitpars, covmat = curve_fit(
			Lorentz,
			freq,
			betrag,
			p0)
			#[f0,fwhm,a0],
			#maxfev=10000
			#)
	print 'fitparameter '+str(fitpars)+'\n\ncovmat '+str(covmat)+'\n\n\n'
	with open('otp_T1.dat','a') as fout:
		fout.write(str(data[7:10])+'\t'+str(np.pi/fitpars[2])+'\n')
	plt.plot(freq,Lorentz(freq,fitpars[0],fitpars[1],fitpars[2]),label=data[7:10])
	p0=fitpars

plt.legend()
plt.show()


