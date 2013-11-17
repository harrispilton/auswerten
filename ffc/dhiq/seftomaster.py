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
import matplotlib.cm as cm
from scipy.optimize import leastsq
from lmfit import  minimize, Parameters, report_errors

def residual(params,xdata,ydata=None,releps=None):
	K_dd=params['K_dd'].value
	tau=params['tau'].value
	beta=params['beta'].value
	if ydata is None:
		return Chi_dd(xdata,K_dd,tau/beta,beta)
	if releps is None:
		return (ydata-Chi_dd(xdata,K_dd,tau/beta,beta))
	return ((ydata-Chi_dd(xdata,K_dd,tau/beta,beta))*releps)
def J_p(nu,tau):
	omega=nu*2*np.pi
	return (tau/(1+(omega*tau)**2))
def J_cd(omega,tau,beta):
	a=omega*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1.+a**2.)**(beta/2.)))
def Chi_dd(omega,K_dd=1.,tau=1.,beta=0.5):
	return omega*K_dd*(J_cd(omega,tau,beta)+4.*J_cd(2.*omega,tau,beta))
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


params= Parameters()
params.add('logK_dd',value=9.0,min=8,max=10)
params.add('K_dd',expr='(10.0**logK_dd)')
params.add('logtau',value=-6.0,min=-12,max=3)
params.add('tau',expr='(10.0**logtau)')

params.add('beta',value=0.9,vary=True,min=0.03,max=1.0)

#text.usetex: True
###
### variablen zuweisen ordner durchfilzen
### 

sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
#c=[]
#for i in np.arange(40):c.append(cm.jet(i/40.))
plt.ion()
plt.figure(1)
ax=plt.axes([0.1,0.1,0.85,0.85])
title=''#raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\nu$ in $MHz$")
plt.ylabel(r"$\chi$ in $s^{-2}$")
plt.xscale('log')
plt.yscale('log')
axcolor = 'lightgoldenrodyellow'

markers=get_markers()
colors=get_colors()

sefdata=[]
temps=[]
brlxs=[]
omegas=[]
chis=[]
r1s=[]
omegas=[]
percerrs=[]
chinorms=[]
omegataus=[]
taus=[]##liste mit log10(tau_strukturrelaxation)
kdds=[]
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	brlx=[]
	omega=[]
	chi=[]
	r1=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
		if liste[0]=='#' or float(liste[2])<0.003:
			pass
		else:
			brlx.append(float(liste[0])*1.e6)
			omega.append(float(liste[0])*1.e6*2.*np.pi)
			r1.append(float(liste[2]))
			chi.append(float(liste[2])*brlx[-1])
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
	brlxs.append(brlx)
	omegas.append(omega)
	chis.append(chi)
	chinorms.append(chi)
	omegataus.append(omega)
	taus.append(0.0)
	r1s.append(r1)
	kdds.append(1e9)
	percerrs.append(percerr)
	#print repr(temp)
	ax.plot(brlx,chi,label=temp+' K',marker=markers.next(),linestyle='None',color=colors.next())

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
for (temp,r1,brlx) in zip(temps,r1s,brlxs):
	with open('r1/ra'+str(temp)+'.dat','w') as fout:
		fout.write('nu(mhz) '+str(temp)+'K\n\n')
		for (r,b) in zip(r1,brlx):
			fout.write(str(b/1.e6)+' '+str(r)+'\n')
for (temp,chi,brlx) in zip(temps,chis,brlxs):
	with open('chi/chi'+str(temp)+'.dat','w') as fout:
		fout.write('nu(mhz) '+str(temp)+'K\n\n')
		for (c,b) in zip(chi,brlx):
			fout.write(str(b/1.e6)+' '+str(c)+'\n')
ax.legend()
plt.draw()
####Normierung der Hoehe:
####aus datensaetzen mit peak im frequenzfenster kann die
####kopplungskonstante bestimmt werden. datensaetze muessen
####vom benutzer gewaehlt werden.

seti=[]
while True:
	seti.append(raw_input("waehle ein datenset mit maximum: "))
	try: seti[-1]=int(seti[-1])
	except ValueError: print 'n zum beenden'
	if seti[-1]=='n': 
		seti.pop()
		break
ks=[]
for i in seti:
	brlxs[i]=np.array(brlxs[i])
	k,tau,beta=1e-5,2e-6,0.9
	out = minimize(residual, params,args=(brlxs[i],chis[i],percerrs[i]),method=('leastsq'))
	result=brlxs[i]+out.residual
	fit = residual(params,brlxs[i])
	ks.append(params['K_dd'].value)
	print params['beta'].value
	print params['tau'].value
####Uebereinanderschieben der Daten:
####die taus koennen von Hand eingegben werden,
####die Strukturrelaxationszeiten werden logarithmisiert in
####einer Liste abgelegt
with open('tau.dat','r') as tauin:
	lines=tauin.readlines()
	if lines.__len__()==taus.__len__():
		beta=params['beta'].value
		for i in range(0,lines.__len__()):
			liste=lines[i].split()
			taus[i]=float(liste[1])
			ax.lines[i].set_xdata([brlx*10**taus[i] for brlx in brlxs[i]])
		plt.autoscale()
		plt.draw()
	else: a=raw_input('laenge der tau stimmt nicht...')


omegatau=np.logspace(-6,4,100)
ax.plot(omegatau,Chi_dd(omegatau))

while True:
	sel=raw_input("Waehle den datenset: ")
	try: 
		int(sel)
		print "aktuelles kdd: "+str(kdds[int(sel)])
		while True:
			kdd=raw_input("neues kdd: ")
			if kdd=='n':break
			beta=params['beta'].value
			kdds[int(sel)]=float(kdd)
			chinorms[int(sel)]=[c/float(kdd) for c in chis[int(sel)]]
			tau=(10**float(taus[int(sel)]))
			ax.lines[int(sel)].set_ydata(chinorms[int(sel)])
			plt.draw()
	except ValueError: print 'n zum beenden'
	if sel=='n':
		break
	if sel=='b':
		beta=float(raw_input('neues beta: '))
		params['beta'].value=beta
		ax.lines[taus.__len__()].set_ydata(Chi_dd(omegatau,1.,1./beta,beta))
		plt.draw()
	if sel=='k':
		k=float(raw_input('neues k: '))
		params['K_dd'].value=k
		beta=params['beta'].value
		for i in range(0,taus.__len__()):
			chinorms[i]=[c/k for c in chis[i]]
			omegataus[i]=[om*10**taus[i] for om in omegas[i]]
			ax.lines[i].set_xdata(omegataus[i])
			ax.lines[i].set_ydata(chinorms[i])
		fit=Chi_dd(omegatau,1.,1./beta,beta)
		ax.lines[taus.__len__()].set_ydata(fit)
		plt.autoscale()
		plt.draw()
		
	if sel=='a':
		minsel=raw_input('schiebe mehrere datensets\n waehle set min:')
		maxsel=raw_input('waehle set max: ')
		logtau=raw_input('um wieviele dekaden sollen die daten geschoben werden? ')
		beta=params['beta'].value
		for i in range(int(minsel),int(maxsel)):
			taus[i]=taus[i]+float(logtau)
			ax.lines[i].set_xdata([om*10**taus[i] for om in omegas[i]])
		ax.autoscale()
		plt.draw()
with open('kdd.dat','w') as fout:
	for i  in range(0,taus.__len__()):
		fout.write(str(temps[i])+' '+str(kdds[i])+'\n')

for i in range(0,taus.__len__()):
	with open('master/ma'+str(temps[i])+' K.dat','w') as fout:
		fout.write('omegatau '+str(temps[i])+'\n\n')
		for (om,ch) in zip(omegataus[i],chinorms[i]):
			fout.write(str(om)+' '+str(ch)+'\n')
with open('master/fit.dat','w') as fout:
	fout.write('beta('+str(params['beta'].value)+') K_dd('+str(params['K_dd'].value)+'\n\n')
	for (om,ch) in zip(omegatau,fit):
		fout.write(str(om)+' '+str(ch)+'\n')


ax.xlabel=(r'$\omega \tau$')
ax.ylabel=(r'$\frac{\chi}{K_(dd)}$')
ax.autoscale()
plt.legend()
plt.draw()


mk=raw_input('ende')
