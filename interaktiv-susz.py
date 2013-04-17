#!/usr/bin/python
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
br=[]
rf=[]
ra1=[]
ch=[]
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
def update(val):
	fin=open(sef[int(picker.val)],'r')
	sefdata=fin.readlines()
	for i in range(0,4):sefdata.pop(0)
	ch=[]
	br=[]
	ra1=[]
	rf=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		br.append(liste[0])
		br=map(float,br)
		ra1.append(liste[2])
		ra1=map(float,ra1)
		rf.append(liste[6])
	print rf[0]	
	slide=10.0**val
#	for i,b in enumerate(br):
#		br[i]=br[i]*10e6*10
#		ch.append(ra1[i]*br[i])
	for i,b in enumerate(br): 
		br[i]=br[i]*10e6*slide
		ch.append(ra1[i]*br[i]/slide)
	plt.plot(br,ch,label=rf[0])
	plt.draw()
def pick(val):
	return val
ax = plt.axes([0.1,0.2,0.55,0.8])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('omega')
plt.ylabel('Chi')

axcolor = 'lightgoldenrodyellow'
axtau_c=plt.axes([0.6,0.1,0.3,0.02],axisbg=axcolor)
stau_c=Slider(axtau_c,'log (tau_c)',-9,9,valinit=np.log10(tau_c))
axpicker=plt.axes([0.1,0.1,0.25,0.02],axisbg=axcolor)
picker=Slider(axpicker,'pick set',0,sef.__len__()-0.01,valinit=0)
sefdata=[]
plt.axes([0.1,0.2,0.55,0.8])

stau_c.on_changed(update)
print picker.on_changed(pick)

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
	plt.plot(brlx,chi,label=relativefile[0])
plt.plot(omega, Chi(omega,1e-6),label='chi mit tau_c =1')
print ax.lines[0]
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
resetax =plt.axes([0.8,0.025,0.1,0.04])
button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')
#set1ax=plt.axes([0.7,0.025,0.1,0.04])
#set1 = Button(set1ax,'set1 relativefile?? somehow',color=axcolor) 
def reset(event):
	stau_c.reset()
button.on_clicked(reset)

plt.axes([0.1,0.2,0.55,0.8])
plt.show()

#vor dem plt.show kann man noch diese 5 zeilen pasten
#rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)

