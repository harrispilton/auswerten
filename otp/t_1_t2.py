#!/usr/bin/python
import glob
import re
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

def lorentz(f,f0,fwhmi,a0):
	return (a0 * (0.5 * fwhm)/((f-f0)**2 + 0.25 * fwhm**2))

files=glob.glob('otp/1d/'+'*.dat')
files.sort()
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
	plt.plot(freq,betrag,label=data,marker='o',linestyle='None')
plt.show()

