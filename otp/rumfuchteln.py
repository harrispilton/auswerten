import numpy as np
from scipy.optimize import curve_fit
def func(x, a, b, c):
	return (a-x**2 + b*x + c)
	# Data (including uncertainties)
x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
y = np.array([-1.78, 4.09, 8.85, 17.9, 26.1, 35.2])
yerr = np.array([0.460, 0.658, 0.528, 1.34, 1.09, 0.786])
p0= [2e5, 3e4,1.0]
def Lorentz(x,f0,fwhm,a0):
	return (a0 * (0.5 * fwhm)/((x * f0)**2.0 + 0.25 * fwhm**2))


f=open('otp/1d/282.5_K.dat','r')
lines=f.readlines()
for i in range(0,2):lines.pop(0)
freq=[]
betrag=[]
for line in lines: 
	liste=line.split()
	freq.append(np.float(liste[0]))
	betrag.append(((np.float(liste[1]))**2+(np.float(liste[2])**2)**0.5))
freq=np.array(freq)
betrag=np.array(betrag)
print type(freq), type(betrag), type(p0),type(p0[0])
# initial guess at parameters
popt, pcov= curve_fit(Lorentz, freq, betrag, p0)
print 'optimal parameters: ', popt
print 'uncertainties of parameters: ', pcov

