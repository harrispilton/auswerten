#!/usr/bin/python
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
beta=0.9
omega=np.logspace(5.0,8.5,20,10)
tau_c=1
K_DD=1.0
delta_sigma_CSA=0.226
def K_CSA(omega):
	a=1.0#2.0/18*(omega * delta_sigma_CSA)**2
	return a
def J(omega,tau_c):
	a=tau_c * 1.0/(1.0+(omega * tau_c)**2)**beta
	return a
def R_1(omega,tau_c):
	a=1/(K_DD * J(omega,tau_c) + K_CSA(omega) * J(omega,tau_c))
	return a

ax = plt.subplot(111)
plt.subplots_adjust(bottom=0.25)
alle=[]
temp=glob.glob('*.dat')
temp.remove('all.dat')
temp.sort()
brlxs=[]
chis=[]
r1=[]
for i, val in enumerate(temp):
	x, y1, y2 = np.loadtxt(temp[i], unpack=True, usecols=[0,1,2])
	brlx = list(x)
	t1 = list(y1)
	r1 = list(y2)
	cbrlx=[]
	wurzelcbrlx=[]
	ct1=[]
	cr1=[]
	chi=[]
	#aussortieren von messpunkten mit gemessenem t1 kleiner als switching time, daher neue listen (damit die laenge passt)
	for ii, val in enumerate(brlx):
		if t1[ii]>0.002:
			cbrlx.append(brlx[ii])
			ct1.append(t1[ii])
			cr1.append(r1[ii])
			chi.append(0)
	for ii, val in enumerate(cbrlx):
		cbrlx[ii]=10e6*cbrlx[ii]
		chi[ii]=10e6*cbrlx[ii]/ct1[ii]	
		wurzelcbrlx.append(cbrlx[ii]**0.5)
	brlxs.append(cbrlx)
	chis.append(chi)
	print wurzelcbrlx
	print cr1
	plt.plot(wurzelcbrlx, cr1,label=temp[i])
#l1, = ax.plot(omega, R_1(omega,tau_c),'r--')
#plt.legend(loc='upper left')
plt.legend()
#plt.yscale('log')
#plt.xscale('log')
#plt.show()

axcolor = 'lightgoldenrodyellow'
axtau_c=plt.axes([0.25,0.1,0.65,0.03],axisbg=axcolor)
stau_c=Slider(axtau_c,'log (tau_c)',-3,6,valinit=np.log10(tau_c))

def update(val):
	l1.set_ydata(R_1(omega,10.0**stau_c.val))
	plt.draw()
stau_c.on_changed(update)
resetax =plt.axes([0.8,0.025,0.1,0.04])
button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')

def reset(event):
	stau_c.reset()
button.on_clicked(reset)
plt.show()
#vor dem plt.show kann man noch diese 5 zeilen pasten
#rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)


#for i,val in enumerate(temp):
#	f = open(temp[i],'r')
#	lines=f.readlines()
#	alle.append[lines],
