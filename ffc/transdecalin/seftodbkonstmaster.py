#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.optimize import brentq 
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.optimize import leastsq
from lmfit import  minimize, Parameters, report_errors, fit_report
from matplotlib.widgets import Slider, Button, CheckButtons

def calc_B():
	mu_0 =1.2566e-6
	h_quer = 6.626e-34/(2.*np.pi)
	gamma_H=2.675e8#/2/np.pi
	N_a=6.022e23
	n_H=18.0
	rho=rho_decalin=.896*1e6#daten fuer cis trans mischung
	M=M_decalin=138.25
	N=n_H*N_a*rho/M
	return np.pi/30.*(1.+4.*(2.**0.5))*(mu_0/4./np.pi * h_quer * gamma_H **2.)**2. * N
def calc_taures(d,r0):
	b=calc_B()
	return (b/(d**1.5*r0))**2.
def calc_konst(d,r0):
	b=calc_B()
	steigung=calc_B()*(d**(-1.5))
	return r0/(steigung**(2./3.))
def calc_r0konst(d,konst):
	steigung=calc_B()*(d**(-1.5))
	return konst*steigung**(2./3.)	
def R_1(sqrtom,R1_0,D):
	B=calc_B()
	return R1_0-(B/(D**1.5))*sqrtom
def J_p(omega,tau):
	return (tau/(1.+(omega*tau)**2.))
def J_cd(omega,tau,beta):
	a=omega*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1.+a**2)**(beta/2.)))
#def Chi_dd(omega,K_dd=1,tau=1,beta=0.5):
#	return omega*3*K_dd*J_cd(omega,tau,beta)
def difvft(T,D0,T0,B):
	return D0*np.exp(B/(T-T0))
def vftres(params,xdata,ydata=None):
	B=vftparams['B'].value
	T0=vftparams['T0'].value
	d=vftparams['D0'].value
	if ydata==None: 
		return difvft(xdata,d,T0,B)
	return (ydata-difvft(xdata,d,T0,B))
def residuals(params,xdata,ydata=None):
	D=params['D'].value
	r0=params['r0'].value
	if ydata==None:
		return R_1(xdata,r0,D)
	return (ydata-R_1(xdata,r0,D))
def update(val):
	konst=skonst.val
	plt.figure(1)
	for i in range(0,temps.__len__()):
		diffs[i]=sdiffs[i].val
		r0s[i]=calc_r0konst(10**diffs[i],konst)
		omtaures=[]
		r1norm=[]
		taures=calc_taures(10**diffs[i],r0s[i])
		for (ome, r) in zip(omegas[i],r1s[i]):
			omtaures.append((ome*taures)**0.5)
			r1norm.append(r/r0s[i])
		ax.lines[i].set_xdata(omtaures)
		ax.lines[i].set_ydata(r1norm)
		konsts.append(calc_konst(10**diffs[i],r0s[i]))
		sqrtoms[i]=omtaures
		r1norms[i]=r1norm
	plt.draw()
def get_colors():
	return itertools.cycle(['g','b','k','c','r','m','0.6'])
def get_markers():
	markers=[]
	for m in plt.Line2D.markers:
		try:
			if len(m)==1 and m !=' ' and m !='|' and m!='_' and m!='x' and m!='.' and m!=',':
				markers.append(m)
		except TypeError:
			pass
	return itertools.cycle(markers)

omega=np.logspace(-3,1.5,200,10)
#K_dd=1e-9
#beta=0.4
#tau_alpha=1
vftparams=Parameters()
vftparams.add('B',value=50,min=10,max=100)
vftparams.add('D0',value=-10,min=-17,max=-5)#expr='10**logD0')
vftparams.add('T0',value=50,min=-80,max=150)
params= Parameters()
params.add('logD',value=-10.0,min=-13.8,max=-9.1)
params.add('D',expr='(10.0**logD)')
params.add('logr0',value=2.,min=-2.5,max=4.)
params.add('r0',expr='(10.0**logr0)')
#params.add('logK_dd',value=8.0,min=7,max=10)
#params.add('K_dd',expr='(10.0**logK_dd)')
#params.add('logtau',value=-6.0,min=-12,max=3)
#params.add('tau',expr='(10.0**logtau)')
#params.add('beta',value=0.5,vary=False,min=0.1,max=1.0)


#plt.figure(2)
#ax=plt.axes([0.1,0.15,0.8,0.8])
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.plot(omega,Chi_dd(omega))
#plt.draw()
#text.usetex: True
###
### variablen zuweisen ordner durchfilzen
### 

sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.85,0.85])
title='Diffusion und so'#raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\sqrt{\omega \tau_{res}}$")
plt.ylabel(r"$R_1$ $[s^{-1}]$")
#plt.xscale('log')
plt.yscale('log')
axcolor = 'lightgoldenrodyellow'
markers=get_markers()#itertools.cycle(['o','s','v','x'])
colors=get_colors()#itertools.cycle(['g','b','k','c','r','m','0.6'])



#plt.figure(2)
#errax=plt.axes([0.1,0.1,0.8,0.8])

sefdata=[]
temps=[]
brlxs=[]
sqrtoms=[]
r1s=[]
r0s=[]
chis=[]
omegas=[]
taus=[]##liste mit log10(tau_strukturrelaxation
diffs=[]
konsts=[]
r1norms=[]
for filename in sef:
	print 'filename: '+filename
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	acolor=colors.next()
	amarker=markers.next()
	omega=[]
	sqrtom=[]
	r1=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
		if liste[0]=='#' or float(liste[1])<0.003:
			pass
		else:
			omega.append(float(liste[0])*1.e6*2.*np.pi)
			sqrtom.append((float(liste[0])*1.e6*2.*np.pi)**0.5)
			r1.append(float(liste[2]))
			percerr.append(float(liste[3]))
			zone.append(int(liste[5]))
			relativefile.append(liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	#wenn einzelne files fehlerhaft sind kann man es einfach durch auskommentieren der folgenden zeilen sehen
	#print 'zone: '+str(zone)
#print 'relativefile: '+relativefile[1]
#print 'zone[sef.index(filename)]' + str(zone[1])
#print 'ZONE=\t'+str(zone[1])+'\n\n'
	temp=sdfdata[sdfdata.index('ZONE=\t'+str(zone[1])+'\r\n')+7]
	temp=temp[6:]
	temp=temp.rstrip()
	temps.append(float(temp))
	ds=[]
	derrs=[]
	iis=[]
	for i in range(r1.__len__()-3,-1,-1):
		minimize(residuals,params,args=(np.array(sqrtom[i:r1.__len__()]),np.array(r1[i:r1.__len__()])))
		iis.append(i)
		ds.append(params['logD'].value)
		derrs.append(params['logD'].stderr)
	
	minni=-1
	minnval=1.e90
	for i in range(1,derrs.__len__()):
		if derrs[i]<minnval:
			minni=iis[i]
			minnval=derrs[i]
	#errax.plot(iis,derrs,label=temp,marker=amarker,color=acolor)
	plt.ylim([0,1.2])
	plt.legend()
	if i == 1:
		pass
	else:
		minimize(residuals,params,args=(np.array(sqrtom[minni:sqrtom.__len__()]),np.array(r1[minni:sqrtom.__len__()])))
		diffs.append(params['logD'].value)
		r0s.append(params['r0'].value)
		fit=residuals(params,np.array(sorted(sqrtom)))
		#print repr(temp)
		omtaures=[]
		r1norm=[]
		konsts.append(calc_konst(10**diffs[-1],r0s[-1]))
		taures=calc_taures(10**diffs[-1],r0s[-1])
		for (ome, r) in zip(omega,r1):
			omtaures.append((ome*taures)**0.5)
			r1norm.append(r/r0s[-1])
		ax.plot(omtaures,r1norm,label=temp+' K',marker=amarker,ms=4.0,color=acolor,linestyle='None')
	#out=minimize(residuals, params,args=(np.array(sqrtom),np.array(r1)))
	omegas.append(omega)
	r1s.append(r1)
	sqrtoms.append(omtaures)
	r1norms.append(r1norm)
	taus.append(0.0)
plt.figure(4)
axdiffs=[]
sdiffs=[]
for i in range(0,temps.__len__()):
	axdiff=plt.axes([0.1,i/(temps.__len__()+2.),0.8,0.02],axisbg=axcolor)
	sdiff=Slider(axdiff,'D('+str(temps[i])+'K)',diffs[i]-2.5,diffs[i]+2,valinit=diffs[i])
	axdiffs.append(axdiff)
	sdiffs.append(sdiff)
#	sdiffs.append(axdiff,
#axdiff=plt.axes([0.1,0.1,0.8,0.05],axisbg=axcolor)
#sdiff=Slider(axdiff,'log D',-15,-9,valinit=-10)
axkonst=plt.axes([0.1,temps.__len__()/(temps.__len__()+2.),0.8,0.05],axisbg=axcolor)
skonst=Slider(axkonst,'log konst',konsts[0]/10.,konsts[0]*10.,valinit=konsts[0])
#axpicker=plt.axes([0.1,(temps.__len__()+1)/(temps.__len__()+2.),0.8,0.05],axisbg=axcolor)
#spicker=Slider(axpicker,'pick set',0,temps.__len__()-0.01,valinit=0)
#spicker.on_changed(pick)
for sdiff in sdiffs:
	sdiff.on_changed(update)
skonst.on_changed(update)
#sdiff.on_changed(update)

with open('D.dat','w') as fout:
	for temp,d in zip(temps,diffs):
		fout.write(str(temp)+' '+str(d)+'\n')

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
ax.legend()
omtau=np.linspace(0.0001,0.9999,105)
ax.plot(omtau**0.5,1.-omtau**0.5,linestyle='--',color='k',label='Modell')
plt.draw()
i=raw_input('next level')
plt.figure(2)
kax=plt.axes([0.1,0.1,0.8,0.8])
konst=(10**diffs[-1]*r0s[-1]+10**diffs[-2]*r0s[-2]+10**diffs[-4]*r0s[-4])/3.
plt.figure(1)
#with open('Jmaster_parameter.dat','r') as fin:
#	lines=fin.readlines()
#	diffs=[]
#	r0s=[]
#	konsts=[]
#	for line in lines:
#		liste=line.split()
#		diffs.append(float(liste[1]))
#		r0s.append(float(liste[2]))
#	print konst
#	r0s=[]
#	for d in diffs:
#		r0s.append(konst/10.**d)
#	for i in range(0,diffs.__len__()):
#		omtaures=[]
#		r1norm=[]
#		taures=calc_taures(10**diffs[i],r0s[i])
#		for (ome, r) in zip(omegas[i],r1s[i]):
#			omtaures.append((ome*taures)**0.5)
#			r1norm.append(r/r0s[i])
#		ax.lines[i].set_xdata(omtaures)
#		ax.lines[i].set_ydata(r1norm)
#		konsts.append(calc_konst(10**diffs[i],r0s[i]))
#		plt.draw()
#		sqrtoms[i]=omtaures
#		r1norms[i]=r1norm
plt.figure(2)
kax.plot(1000./np.array(temps),konsts)
plt.figure(3)
dax=plt.axes([0.1,0.1,0.8,0.8])
dax.scatter(temps,diffs)
minimize(vftres,vftparams,args=(np.array(temps),np.array(diffs)))
fit=vftres(vftparams,np.array(temps))
print report_errors(vftparams)
dax.plot(temps,fit)
plt.draw()
while True:
	sel=raw_input('waehle set:  ')
	if sel=='n':break
	sel=int(sel)
	print 'aktuelles r1: '+str(r0s[sel])+'\naktuelles d:  '+str(diffs[sel])
	rod=raw_input('r oder d')
	if rod =='k':
		kneu=float(raw_input('neue konst: '))
		for i in range(0,diffs.__len__()):
			konsts[i]=kneu
			konst=kneu
			r0s[i]=kneu/(10.**diffs[i])
			omtaures=[]
			r1norm=[]
			b=calc_B()
			#taures=(b/((10**diffs[sel])**1.5*r0neu))**2.
			taures=calc_taures(10**diffs[i],r0s[i])
			for (ome, r) in zip(omegas[i],r1s[i]):
				omtaures.append((ome*taures)**0.5)
				r1norm.append(r/r0s[i])
			ax.lines[i].set_xdata(omtaures)
			ax.lines[i].set_ydata(r1norm)
			plt.draw()
			sqrtoms[i]=omtaures
			r1norms[i]=r1norm
	if rod =='r':
		r0neu=float(raw_input('neues r0: '))
		omtaures=[]
		r1norm=[]
		taures=calc_taures(10**diffs[sel],r0neu)
		for (ome, r) in zip(omegas[sel],r1s[sel]):
			omtaures.append((ome*taures)**0.5)
			r1norm.append(r/r0neu)
		ax.lines[sel].set_xdata(omtaures)
		ax.lines[sel].set_ydata(r1norm)
		konsts[sel]=calc_konst(diffs[sel],r0neu)
		plt.figure(2)
		plt.cla()
		kax.plot(1000./np.array(temps),konsts)
		plt.draw()
		plt.figure(1)
		plt.draw()
		r0s[sel]=r0neu
		sqrtoms[sel]=omtaures
		r1norms[sel]=r1norm
	elif rod =='d':
		dneu=float(raw_input('neues d: '))
		r0s[sel]=konst/10**dneu
		omtaures=[]
		r1norm=[]
		b=calc_B()
		steigung=calc_B()*((10**dneu)**(-1.5))
		konsts[sel]=(r0s[sel]/(steigung**(2./3.)))
		taures=(b/((10**dneu)**1.5*r0s[sel]))**2.
		for (ome, r) in zip(omegas[sel],r1s[sel]):
			omtaures.append((ome*taures)**0.5)
			r1norm.append(r/r0s[sel])
		ax.lines[sel].set_xdata(omtaures)
		ax.lines[sel].set_ydata(r1norm)
		plt.figure(2)
		plt.cla()
		kax.plot(1000./np.array(temps),konsts)
		plt.draw()
		plt.figure(1)
		plt.draw()
		diffs[sel]=dneu
	else:
		pass
for (ome,r1,temp) in zip(sqrtoms,r1norms,temps):
	with open('dbkonstma/dbkonst'+str(temp)+'.dat','w') as fout:
		fout.write(str(temp)+' '+str(temp)+'\n\n')
		for (o,r) in zip(ome,r1):
			fout.write(str(o)+' '+str(r)+'\n')

with open('Jmaster_parameter.dat','w') as fout:
		for (t,d,r,k) in zip(temps,diffs,r0s,konsts):
			fout.write(str(t)+' '+str(d)+' '+str(r)+' '+str(calc_B())+' '+str(k)+'\n')


oksdfj=raw_input('ente')
