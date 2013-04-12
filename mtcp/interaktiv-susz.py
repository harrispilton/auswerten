#!/usr/bin/python
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
beta=0.9
omega=np.logspace(5.0,8.5,20,10)
tau_c=1.0e-6
K_DD=1.0e9
delta_sigma_CSA=0.226
sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
def J(omega,tau_c):
	return tau_c/(1. + omega **2 * tau_c **2)**beta
def Chi(omega,tau_c):
	return K_DD*omega*J(omega,tau_c)
def Transform(a, brlx=[],r1=[]):
	print r1
	print type(r1)
	for rate in r1: r1[i]=rate*a*brlx[1]
	return r1
plt.axes([0.1,0.1,0.55,0.8])
sefdata=[]
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	chi=[]
	brlx=[]
	t1=[]
	r1=[]
	percerr=[]
	abserr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		brlx.append(liste[0])
		brlx=map(float,brlx)
		t1.append(liste[1])
		t1=map(float,t1)
		r1.append(liste[2])
		r1=map(float,r1)
		percerr.append(liste[3])
		percerr=map(float,percerr)
		abserr.append(liste[4])
		abserr=map(float,abserr)
		zone.append(liste[5])
		zone=map(int,zone)
		relativefile.append(liste[6])
	for i,b in enumerate(brlx):
		brlx[i]=brlx[i]*10e6
		chi.append(r1[i]*brlx[i])
	print type(r1)	
	if filename == "mtcp_200K.sef": 
		print r1[1]
		Transform(brlx,chi,2.0)
	plt.plot(brlx,chi,label=relativefile[0])
print Chi(20e7,1e-9)
plt.plot(omega, Chi(omega,1e-6),label='chi mit tau_c =1')

plt.legend(loc='center left',bbox_to_anchor=(1,0.5))

#
#axcolor = 'lightgoldenrodyellow'
#axtau_c=plt.axes([0.25,0.1,0.3,0.02],axisbg=axcolor)
#
#stau_c=Slider(axtau_c,'log (tau_c)',-3,6,valinit=np.log10(tau_c))
#
#def update(val):
#	l1.set_ydata(Chi(omega,10.0**stau_c.val))
#	plt.draw()
#stau_c.on_changed(update)
#resetax =plt.axes([0.8,0.025,0.1,0.04])
#button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')
#
#def reset(event):
#	stau_c.reset()
#button.on_clicked(reset)
#



#vor dem plt.show kann man noch diese 5 zeilen pasten
#rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)


plt.xscale('log')
plt.yscale('log')
plt.xlabel('omega')
plt.ylabel('Chi')
plt.show()
