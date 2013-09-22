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

plt.ion()
plt.figure(1)
ax=plt.axes([0.15,0.15,0.8,0.8])
#ax.set_xscale='lin'
#ax.set_yscale='log'
plt.yscale('log')
satrec=open('satrec.dat','r')
lines=satrec.readlines()
tempsatrec=[]
t1satrec=[]
for line in lines:
	liste = line.split()
	tempsatrec.append(np.float(liste[0]))
	t1satrec.append(np.float(liste[1]))
invrec=open('invrec.dat','r')
lines=invrec.readlines()
tempinvrec=[]
t1invrec=[]
for line in lines:
	liste = line.split()
	tempinvrec.append(np.float(liste[0]))
	t1invrec.append(np.float(liste[1]))

t2=open('otp_T2.dat','r')
lines=t2.readlines()
tempecho=[]
t2echo=[]
for line in lines:
	liste = line.split()
	tempecho.append(np.float(liste[0]))
	t2echo.append(np.float(liste[1]))
ax.plot(tempecho,t2echo,label='T2 aus Solid Echo', linestyle='None',marker='o')
ax.plot(tempinvrec,t1invrec,label='T1 aus Inversion Recovery',linestyle='None',marker='o')
ax.plot(tempsatrec,t1satrec,label='T1 aus Saturation Recovery',linestyle='None',marker='o')
plt.legend()
plt.draw()
plt.show()
i=raw_input('ende')
