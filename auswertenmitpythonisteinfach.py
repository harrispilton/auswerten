#!/usr/bin/python
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
beta=0.6
omega=np.logspace(5.0,8.5,20,10)
tau_c=1e30
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
for i, val in enumerate(temp):
	x, y1, y2 = np.loadtxt(temp[i], unpack=True, usecols=[0,1,2])
	brlx = list(x)
	t1 = list(y1)
	r1 = list(y2)
	cbrlx=[]
	ct1=[]
	cr1=[]
	chi=[]
	for ii, val in enumerate(brlx):
		if t1[ii]>0.002:
			cbrlx.append(brlx[ii])
			ct1.append(t1[ii])
			cr1.append(r1[ii])
			chi.append(0)
	for ii, val in enumerate(cbrlx):
		cbrlx[ii]=10e6*cbrlx[ii]
		chi[ii]=10e6*cbrlx[ii]/ct1[ii]	
	brlxs.append(cbrlx)
	chis.append(chi)
	plt.plot(cbrlx, chi,label=temp[i])

l1, = ax.plot(omega, R_1(omega,tau_c),'r--')
plt.legend(loc='upper left')
plt.yscale('log')
plt.xscale('log')
#plt.show()

axcolor = 'lightgoldenrodyellow'
axtau_c=plt.axes([0.25,0.1,0.65,0.03],axisbg=axcolor)
stau_c=Slider(axtau_c,'log (tau_c)',29,50,valinit=np.log10(tau_c))

def update(val):
	l1.set_ydata(R_1(omega,10.0**stau_c.val))
	plt.draw()
stau_c.on_changed(update)
resetax =plt.axes([0.8,0.025,0.1,0.04])
button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')

def reset(event):
	stau_c.reset()
button.on_clicked(reset)
rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)
plt.show()


#axtau_c = plt.axes([0.25,0.1,0.65,0.03])
#stau_c = Slider(axtau_c, 'tau_c', 0.00001, 10.0, valinit=tau_c)
#button = Button(plt.axes([0.8,0.025,0.1,0.04]), 'Reset',hovercolor='0.6')
#def update(val):
#    tau_c = samp.val
#    l.set_ydata(R_1(omega,tau_c))
#    plt.draw()
#stau_c.on_changed(update)
#
#resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
#def reset(event):
#    stau_c.reset()
#button.on_clicked(reset)
#rax = plt.axes([0.025, 0.5, 0.15, 0.15])
#def colorfunc(label):
#    l.set_color(label)
#    plt.draw()
#
#
#for i,val in enumerate(temp):
#	f = open(temp[i],'r')
#	lines=f.readlines()
#	alle.append[lines],
