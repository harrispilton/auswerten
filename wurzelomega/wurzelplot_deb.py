#!/usr/bin/python
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
ax = plt.subplot(111)
x1, y1, y2 = np.loadtxt('dhiq 280K.dat', unpack=True, usecols=[0,1,2])
brlx = list(x1)
t1 = list(y1)
cr1= list(y2)
wurzelcbrlx=brlx
for i in range(brlx.__len__()):	wurzelcbrlx[i] = brlx[i]**0.5
for i in range(8): 
	del wurzelcbrlx[-1]
	print wurzelcbrlx
	del cr1[-1]
#print wurzelcbrlx
#print cr1
def func(x,m,t):
	return t-x*m
p0=[1.5,7.0]
popt = np.polyfit(wurzelcbrlx, cr1,1)
print popt[0], popt[1]
plt.plot(wurzelcbrlx,cr1, 'ko')
x=np.linspace(0,4.5,100)
plt.plot(x, popt[0]*x+popt[1], 'r-')
plt.show()
#fitfunc = lambda p, x: p[0]*x+p[1]
#errfunc = lambda p, x, y: fitfunc(p,x)-y
#p0=[-1.0,1.0]
#p1, success =optimize.leastsq(errfunc,p0[:], args=[wurzelcbrlx,cr1])

