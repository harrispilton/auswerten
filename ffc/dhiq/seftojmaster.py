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

def calc_B():
	mu_0 =1.2566e-6
	h_quer = 6.626e-34/(2.*np.pi)
	gamma_H=2.675e8#/2/np.pi
	N_a=6.022e23
	n_H=21.0
	rho=rho_mTCP=1.15*1e6
	M=M_mTCP=368.4
	N=n_H*N_a*rho/M
	return np.pi/30.*(1.+4.*(2.**0.5))*(mu_0/4./np.pi * h_quer * gamma_H **2.)**2. * N

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
def residuals(params,xdata,ydata=None):
	D=params['D'].value
	r0=params['r0'].value
	if ydata==None:
		return R_1(xdata,r0,D)
	return (ydata-R_1(xdata,r0,D))
def get_colors():
	return itertools.cycle(['g','b','k','c','r','m','0.6'])
def get_markers():
	return itertools.cycle(['o','s','v','x','^'])

omega=np.logspace(-3,1.5,200,10)
#K_dd=1e-9
#beta=0.4
#tau_alpha=1
params= Parameters()
params.add('logD',value=-9.0,min=-14.5,max=-8.5)
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
title=r'Spektraldichte Masterkurve'#raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\sqrt{\omega \tau_{res}}$")
plt.ylabel(r"$R_1$ $[s^{-1}]$")
#plt.xscale('log')
#plt.yscale('log')
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
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	omega=[]
	sqrtom=[]
	r1=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
		if liste[0]=='#':
			pass
		else:
		#	liste = re.findall(r"[\w.][\f]+",data)
			omega.append(float(liste[0])*1.e6*2.*np.pi)
			sqrtom.append((float(liste[0])*1.e6*2.*np.pi)**0.5)
			r1.append(float(liste[2]))
			percerr.append(float(liste[3]))
			zone.append(int(liste[5]))
			relativefile.append(liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	print 'filename: '+filename
	#wenn einzelne files fehlerhaft sind kann man es einfach durch auskommentieren der folgenden zeilen sehen
	#print 'zone: '+str(zone)
#print 'relativefile: '+relativefile[1]
#print 'zone[sef.index(filename)]' + str(zone[1])
#print 'ZONE=\t'+str(zone[1])+'\n\n'
	temp=sdfdata[
		sdfdata.index(
			'ZONE=\t'+str(zone[
				1])+'\r\n')+7]
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
	minni=1
	minnval=1.e90
	for i in range(1,derrs.__len__()):
		if derrs[i]<minnval:
			minni=iis[i]
			minnval=derrs[i]
	acolor=colors.next()
	amarker=markers.next()
	#errax.plot(iis,derrs,label=temp,marker=amarker,color=acolor)
	plt.ylim([0,1])
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
		b=calc_B()
		taures=(b/(params['D']**1.5*params['r0'].value))**2.
		for (ome, r) in zip(omega,r1):
			omtaures.append((ome*taures)**0.5)
			r1norm.append(r/params['r0'].value)
		ax.plot(omtaures,r1norm,label=temp+' K',marker=amarker,ms=4.0,color=acolor,linestyle='None')
	#out=minimize(residuals, params,args=(np.array(sqrtom),np.array(r1)))
	omegas.append(omega)
	r1s.append(r1)
	sqrtoms.append(sqrtom)
	taus.append(0.0)

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
ax.legend()
omtau=np.linspace(0.01,0.99,100)
ax.plot(omtau**0.5,1.-omtau**0.5,linestyle='--',color='k',label='Modell')
#ax.plot(omtau**0.5,R_1(omtau**0.5,1,calc_B()),linestyle='--',color='k',label='fit')
plt.draw()
oksdfj=raw_input('ente')
