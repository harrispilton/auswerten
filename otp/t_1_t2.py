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
	return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))

markers = itertools.cycle(['o','s','v'])
files=glob.glob('otp/1d/'+'*K.dat')
files.sort()
x=np.linspace(-1e6,1e6,1e5)
p0=[2.5e5,1e3,1]
#a0_a=[]
#fwhm_a=[]
#f0_a=[]
#for freq in fr:
#	a0_a.append(a0)
#	fwhm_a.append(fwhm)
#	f0_a.append(f0)
for data in files:
	f=open(data,'r')
	lines=f.readlines()
	for i in range(0,2): lines.pop(0)
	freq=[]
	betrag=[]
	for line in lines:
		liste=line.split()
		freq.append(liste[0])
		betrag.append(((float(liste[1]))**2+(float(liste[2])**2)**0.5))
	plt.figure(1)
	plt.plot(freq,betrag,label=str(data),marker=markers.next(),linestyle='None')	
	#plt.plot(fr,Lorentz(fr,f0,fwhm,a0))
	freq=np.asarray(freq).ravel()
	betrag=np.asarray(betrag).ravel()
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


plt.legend()
plt.show()


