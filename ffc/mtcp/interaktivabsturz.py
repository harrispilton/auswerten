#!/usr/bin/python
import glob
import re
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
#text.usetex: True
###
### variablen zuweisen ordner durchfilzen
### 
br=[]
rf=[]
ra1=[]
ch=[]
beta=0.9
omega=np.logspace(4.0,7.301,20,10)
tau_c=1.0e-6
K_DD=1.0e9
delta_sigma_CSA=0.226
sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
verschiebefaktoren=[]
verschiebetemperaturen=[]
for i in range(0,sef.__len__()): 
	verschiebetemperaturen.append(0)
	verschiebefaktoren.append(0)
wurzelomega=[]
for om in omega: wurzelomega.append(om**0.5)
###
### funktionen definieren
###

##die spektraldichte J und die Suszibilitaet Chi
def J(omega,tau_c,beta):
	return tau_c/(1. + omega **2 * tau_c **2)**beta
def Chi(omega,tau_c,beta,K_DD):
	return K_DD*omega*J(omega,tau_c,beta)
##die auswertefunktion fuer die Diffusion und die rate. mal sehen...
def R_1(omega,R1_0,D):
	mu_0 =1.2566e-6
	h_quer = 6.626e-34/(2.*np.pi)
	gamma_H=2.675e8#/2/np.pi
	N_a=6.022e23
	n_H=21.0
	rho=rho_mTCP=1.15*1e6
	M=M_mTCP=368.4
	N=n_H*N_a*rho/M
	B=np.pi/30.*(1.+4.*(2.**0.5))*(mu_0/4./np.pi * h_quer * gamma_H **2)**2 * N
	#print N, omega,R1_0,D,B,R1_0-B/(D**1.5) *omega**0.5
	return R1_0-B/(D**1.5) *omega**0.5
##die verschiebefunktion fuer die suszibilitaet
def update(val):
	fin=open(sef[int(picker.val)],'r')
	sefdata=fin.readlines()
	for i in range(0,4):sefdata.pop(0)
	ch=[]
	br=[]
	ra1=[]
	zone=[]
	rf=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		br.append(liste[0])
		br=map(float,br)
		ra1.append(liste[2])
		ra1=map(float,ra1)
		zone.append(liste[5])
		zone=map(int,zone)
		rf.append(liste[6])
	slide=10.0**val
	for i,b in enumerate(br): 
		br[i]=br[i]*1e6
		ch.append(ra1[i]*br[i])
		br[i]=br[i]*slide
	#for line in ax.lines: print line
	ax.lines[int(picker.val)*2+1].set_xdata(br)
	plt.draw()
	verschiebefaktoren[int(picker.val)]=slide
	verschiebetemperaturen[int(picker.val)]=temps[int(picker.val)]
	plt.figure(4)
	plt.cla()
	plt.plot(verschiebetemperaturen,verschiebefaktoren)
	plt.draw()
	return slide
#	fin=open(sef[int(picker.val)],'r')
#	sefdata=fin.readlines()
#	for i in range(0,4):sefdata.pop(0)
#	ch=[]
#	br=[]
#	ra1=[]
#	rf=[]
#	for data in sefdata: 
#		liste=data.split()
#	#	liste = re.findall(r"[\w.][\f]+",data)
#		br.append(liste[0])
#		br=map(float,br)
#		ra1.append(liste[2])
#		ra1=map(float,ra1)
#		rf.append(liste[6])
#	slide=10.0**val
#	for i,b in enumerate(br): 
#		br[i]=br[i]*10e6
#		ch.append(ra1[i]*br[i])
#		br[i]=br[i]*slide
#	ax.lines[int(picker.val)].set_xdata(br)
#	plt.draw()
## den pick gibts nur anstandshalber
def pick(val):
	return val
def r_ref(r):
	br=[]
	ra1=[]
	rf=[]
	d=10**float(sd0.val)
	for i,om in enumerate(omega): 
		br.append(omega[i]**0.5)
		ra1.append(R_1(omega[i],r,d))
	wurzelax.lines[wurzelax.lines.__len__()-1].set_ydata(ra1)
	plt.draw()
def d0(d):
	br=[]
	ra1=[]
	rf=[]
	r=float(sr0.val)
	d=10**d
	print d
	for i,om in enumerate(omega): 
		br.append(omega[i]**0.5)
		ra1.append(R_1(omega[i],r,d))
	wurzelax.lines[wurzelax.lines.__len__()-1].set_ydata(ra1)
	plt.draw()
	return d
def normchi(K):
	plt.figure(1)
	for line in ax.lines:
		linex=line.get_xdata()
		liney=line.get_ydata()
		for y in liney: y=y/(10**K)
		line.set_ydata(liney)
	plt.draw()
def reset(event):
	stau_c.reset()

plt.figure(1)
ax = plt.axes([0.1,0.2,0.55,0.7])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('omega')
plt.ylabel('Chi')

axcolor = 'lightgoldenrodyellow'
axtau_c=plt.axes([0.6,0.1,0.3,0.02],axisbg=axcolor)
stau_c=Slider(axtau_c,'log (tau_c)',-7,2,valinit=np.log10(tau_c))
axpicker=plt.axes([0.1,0.1,0.25,0.02],axisbg=axcolor)
picker=Slider(axpicker,'pick set',0,sef.__len__()-0.01,valinit=0)
resetax =plt.axes([0.8,0.025,0.1,0.04])
button = Button(resetax,'reset',color=axcolor,hovercolor='0.975')
axK = plt.axes([0.1,0.05,0.3,0.02],axisbg=axcolor)
sK=Slider(axK,'K',8,12,valinit=9)

plt.figure(2)
wurzelax=plt.axes([0.1,0.1,0.8,0.8])
axr0=plt.axes([0.05,0.02,0.6,0.02],axisbg=axcolor)
sr0=Slider(axr0,'r0',0.05,4,valinit=1.0)
axD=plt.axes([0.7,0.02,0.2,0.02],axisbg=axcolor)
sd0=Slider(axD,'D',-15,-7,valinit=-11)

plt.figure(3)
##wird spaeter bemalt

plt.figure(4)
fakchiax=plt.axes([0.15,0.15,0.8,0.8])
plt.title('verschiebefaktoren in der suszeptiblitaet')
plt.xlabel('T')
plt.ylabel('schiebefaktoren a.u.')

##interaktion mit der gui
picker.on_changed(pick)
button.on_clicked(reset)
stau_c.on_changed(update)
sK.on_changed(normchi)
sr0.on_changed(r_ref)
sd0.on_changed(d0)

sefdata=[]
temps=[]
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
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	temp=sdfdata[
			sdfdata.index(
				'ZONE=\t'+str(zone[
					sef.index(filename)])+'\r\n')+7]
	temp=temp[6:]
	temp=temp.rstrip()
	temps.append(float(temp))
	#print repr(temp)
	for i,b in enumerate(brlx):
		brlx[i]=brlx[i]*1e6
		chi.append(r1[i]*brlx[i])
	plt.figure(1)
	plt.title(relativefile[0][0:-9])
	plt.axes([0.1,0.2,0.55,0.7])
	#plt.plot(brlx,map(lambda x:Chi(x,1e-6,0.7,1e8),brlx))
	brlx=np.array(brlx)
	chi=np.array(chi)
	fitpars, covmat = curve_fit(
			Chi,
			brlx,
			chi,
			p0=[1e-6,0.7,1e10],
			maxfev=10000)
	print 'fitparamer temperatur (tau,beta,kopplungskonstante)'+temp+str(fitpars)
	plt.plot(brlx,
			map(lambda x:Chi(x,fitpars[0],fitpars[1],fitpars[2]),brlx),
			#label='Chi '+temp+' K'
			)
	plt.plot(brlx,chi,
			label=temp+' K',
			marker='o',linestyle='None')
	plt.figure(2)
	plt.title('Wurzel')
	wurzelax=plt.axes([0.1,0.1,0.8,0.8])
	for i, b in enumerate(brlx):
		brlx[i]=brlx[i]**0.5
	plt.plot(brlx,r1,label=temp+' K')
	plt.figure(3)
	plt.title('Rate')
	plt.xscale('log')
	plt.yscale('log')
	for i,b in enumerate(brlx):
		brlx[i]=brlx[i]**2 #wir hatten die wurzel gezogen
	plt.plot(brlx,r1,label=relativefile[0])

plt.figure(4)

print (map(lambda x: x**0.5,omega),map(lambda y:R_1(y,2,1e-10),wurzelomega))

plt.figure(2)
plt.plot(wurzelomega,map(lambda x: R_1(x,20,10e-9),omega))
plt.legend()
plt.figure(1)
plt.plot(omega, Chi(omega,1e-6,0.7,1e8),label='chi mit tau_c =1e-6')
plt.plot(omega, Chi(omega,1e-8,0.7,1e8),label='chi mit tau_c =1e-8')
plt.plot(omega, Chi(omega,1e-4,0.7,1e8),label='chi mit tau_c =1e-4')
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig('suscibility',dpi=300,orientation='landscape')
plt.figure(2)
plt.savefig('wurzel',dpi=300,orientation='landscape')
plt.figure(3)
plt.savefig('rate',dpi=300,orientation='landscape')
plt.figure(4)
plt.savefig('verschiebeparameter',dpi=300,orientation='landscape')
plt.show()


#set1ax=plt.axes([0.7,0.025,0.1,0.04])
#set1 = Button(set1ax,'set1 relativefile?? somehow',color=axcolor) 

#for i in range(0,ax.lines.__len__()-1): print ax.lines[i]
#vor dem plt.show kann man noch diese 5 zeilen pasten
#rax = plt.axes([0.025,0.5,0.15,0.15],axisbg=axcolor)
#def colorfunc(label):
#	l.set_color(label)
#	plt.draw()
#radio.on_clicked(colorfunc)
#def J(x,tau_c,beta):
#	return tau_c /((1.0 + ( x * tau_c )**2 )** beta)
#plt.figure(4)
#plt.plot(omega,J(omega,1e-8,0.5))
#def Chi(x,tau_c,K_DD,beta):
#	return 1.08 * 1.08 * J(x,tau_c,beta)

