#!/usr/bin/python
#import matplotlib
import glob
#matplotlib.use('WXAgg')
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
#from matplotlib.font_manager import FontProperties
#br=[]
#rf=[]
#ra1=[]
#ch=[]
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
	fin=open(sef[int(spicker.val)],'r')
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
	slide=10.0**val
#	for i,b in enumerate(br):
#		br[i]=br[i]*10e6*10
#		ch.append(ra1[i]*br[i])
	for i,b in enumerate(br): 
		br[i]=br[i]*10e6
		ch.append(ra1[i]*br[i])
		br[i]=br[i]*slide
	#plt.plot(br,ch,label=rf[0])
	ax.lines[int(spicker.val)].set_xdata(br)
	plt.draw()
def pick(val):
	return val
def reset(event):
	stau_c.reset()

ax=plt.axes([0.1,0.1,0.55,0.85])
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('omega')
plt.ylabel('Rate in 1/s')
plt.setp(ax.get_xticklabels())
ax.autoscale_view(True,True,True)
#ax2=plt.axes([0.1,0.2,0.55,0.35])
#plt.xscale('log')
#plt.yscale('log')
##plt.xlabel('omega')
#plt.ylabel('Chi')
#ax2.sharex=ax
#ax2.sharey=ax
#plt.setp(ax2.get_xticklabels(),visible=False)


#axcolor = 'lightgoldenrodyellow'
#
#resetax =plt.axes([0.8,0.025,0.1,0.04])
#button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')
#
#axtau_c=plt.axes([0.6,0.1,0.3,0.02],axisbg=axcolor)
#stau_c=Slider(axtau_c,'log (tau_c)',-7,7,valinit=np.log10(tau_c))
#
#axpicker=plt.axes([0.1,0.1,0.25,0.02],axisbg=axcolor)
#spicker=Slider(axpicker,'pick set',0,sef.__len__()-0.01,valinit=0)
#
#stau_c.on_changed(update)
#spicker.on_changed(pick)
#button.on_clicked(reset)

ax
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
		print i, brlx[i]
		brlx[i]=brlx[i]*1e6
		chi.append(r1[i]*brlx[i])
		print i, brlx[i]
		brlx[i]=brlx[i]**0.5
		print i, brlx[i]
	ax.plot(brlx,r1,label=relativefile[0])
	#ax2.plot(brlx,chi)	

ax.legend(loc='center left',bbox_to_anchor=(1.,0.5))
#set1ax=plt.axes([0.7,0.025,0.1,0.04])
#set1 = Button(set1ax,'set1 relativefile?? somehow',color=axcolor) 

plt.show()

#vor dem plt.show kann man noch diese 5 zeilen pasten
#rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)
