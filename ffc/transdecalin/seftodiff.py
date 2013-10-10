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
	n_H=18.0
	rho=rho_decalin=.896*1e6#daten fuer cis trans mischung
	M=M_decalin=138.25
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
params= Parameters()
params.add('logD',value=-11.2,min=-11.95,max=-8.55)
params.add('D',expr='(10.0**logD)')
#params.add('steigung',expr='(calc_B()/((10.0**logD)**1.5))')
#params.add('steigung',expr='(r0/((10.0**logD)**(2.0/3.0)))')
params.add('logr0',value=2.,min=-1.0,max=4.)
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
plt.xlabel(r"$\sqrt{\omega}$")
plt.ylabel(r"$R_1$ $[s^{-1}]$")
#plt.xscale('log')
#plt.yscale('log')
axcolor = 'lightgoldenrodyellow'

markers=get_markers()#itertools.cycle(['o','s','v','x'])
colors=get_colors()#itertools.cycle(['g','b','k','c','r','m','0.6'])

insetax=plt.axes([0.4,0.7,0.2,0.2])
plt.ylabel(r"$lg(D)$")
plt.xlabel(r"$\frac{1000}{T}$ $[\frac{1}{K}]$")
plt.xscale('linear')
plt.yscale('linear')


plt.figure(2)
errax=plt.axes([0.1,0.1,0.8,0.4])
dmax=plt.axes([0.1,0.5,0.4,0.4])


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
ms=[]
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	acolor=colors.next()
	amarker=markers.next()
	brlx=[]
	sqrtom=[]
	r1=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		if liste[0]=='#' or float(liste[3])>100:
			pass
		else:
			brlx.append(float(liste[0])*1.e6*2.*np.pi)
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
	steigungs=[]
	iis=[]
	for i in range(r1.__len__()-3,-1,-1):
		minimize(residuals,params,args=(np.array(sqrtom[i:r1.__len__()]),np.array(r1[i:r1.__len__()])))
		iis.append(i)
		ds.append(params['logD'].value)
		derrs.append(params['logD'].stderr)
		steigungs.append(calc_B()/(params['D'].value**1.5))
	minni=1
	minnval=1.e90
	plt.figure(2)
	dmax.plot(iis,params['r0'].value/(np.array(steigungs)**(2./3.)),color=acolor,marker=amarker)
	dmax.autoscale()
	plt.xscale('log')
	plt.yscale('log')
	#plt.draw()
	#raw_input('next')
	#dmax.cla()
	for i in range(1,derrs.__len__()):
		if derrs[i]<minnval:
			minni=iis[i]
			minnval=derrs[i]
	plt.figure(2)
	errax.plot(iis,derrs,label=temp,marker=amarker,color=acolor)
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
		plt.figure(1)
	
		ax.plot(sqrtom,r1,label=temp+' K',marker=amarker,ms=4.0,color=acolor,linestyle='None')
		ax.plot(sorted(sqrtom),fit,linestyle='--',color=acolor)
		insetax.plot(1000./(float(temp)),params['logD'].value,marker=amarker,color=acolor)
		#out=minimize(residuals, params,args=(np.array(sqrtom),np.array(r1)))
		brlxs.append(brlx)
	r1s.append(r1)
	sqrtoms.append(sqrtom)
	taus.append(0.0)

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
ax.legend()
with open('D.dat','w') as fout:
	for temp,d in zip(temps,diffs):
		print (str(temp)+' '+str(10**d)+'\n')
		fout.write(str(temp)+' '+str(d)+'\n')
plt.draw()
oksdfj=raw_input('ente')
##while True:
##	selecttofit=raw_input("waehle die datensets mit alphapeak (0-"+str(temps.__len__())+") abbrechen mit n:  ")
##	if selecttofit==n:break
##	brlxs[int(selecttofit)]
#####Normierung der Hoehe:
#####aus datensaetzen mit peak im frequenzfenster kann die
#####kopplungskonstante bestimmt werden. datensaetze muessen
#####vom benutzer gewaehlt werden.
#
##k=raw_input("Schaetze die Kopplungskonstante:  ") 
##try:
#	#k=float(k)
#seti=[]
#while True:
#	seti.append(raw_input("waehle ein datenset mit maximum: "))
#	try: seti[-1]=int(seti[-1])
#	except ValueError: print 'n zum beenden'
#	if seti[-1]=='n': 
#		seti.pop()
#		break
#ks=[]
#for i in seti:
#	brlxs[i]=np.array(brlxs[i])
#	k,tau,beta=1e-5,2e-6,0.5
#	out = minimize(residual, params,args=(brlxs[i],chis[i]),method=('leastsq'))
#	result=brlxs[i]+out.residual
#	fit = residual(params,brlxs[i])
#	print fit
#	print result
#	ks.append(params['K_dd'].value)
#	print params['beta'].value
#	print params['tau'].value
#	
##	report_errors(params)
##	plt.figure(3)
##
##	fitax=plt.axes([0.1,0.1,0.8,0.8])
##	fitax.set_xscale('log')
##	fitax.set_yscale('log')
##	plt.plot(brlxs[i],Chi_dd(brlxs[i],params['K_dd'].value,params['tau'].value,0.5))
##	plt.plot(brlxs[i],chis[i])
##	plt.show()
#
##	mk=raw_input('fuer naechster fit ente druecken')
#K_dd=0.0
#for k in ks:
#	K_dd=K_dd+k
#K_dd=K_dd/ks.__len__()
#for i in range(0,chis.__len__()):
#	chis[i]=[chi/K_dd for chi in chis[i]]
#
#####Uebereinanderschieben der Daten:
#####die taus koennen von Hand eingegben werden,
#####Fits bringen hier niemanden weiter
#####Im Folgenden kann man einzelne Daten per Eingabedialog verschieben
#####die Strukturrelaxationszeiten werden logarithmisiert in
#####einer Liste abgelegt
#while True:
#	tauin=open('tau.dat','r')
#	lines=tauin.readlines()
#	if lines.__len__()==taus.__len__():
#		for i in range(0,lines.__len__()):
#			liste=lines[i].split()
#			taus[i]=float(liste[1])
#			ax.lines[i].set_xdata([brlx*10**taus[i] for brlx in brlxs[i]])
#		plt.draw()
#	else: a=raw_input('laenge der tau stimmt nicht...')
#	break
#
#while True:
#	sel=raw_input("Waehle den datenset: ")
#	try: 
#		int(sel)
#		print "aktuelles log tau: "+str(taus[int(sel)])
#		while True:
#			logtau=raw_input("neues tau: ")
#			if logtau=='n':break
#			tau=10**float(logtau)
#			taus[int(sel)]=float(logtau)
#			ax.lines[int(sel)].set_xdata([brlx*tau for brlx in brlxs[int(sel)]])
#			plt.draw()
#	except ValueError: print 'n zum beenden'
#	if sel=='n':break
#	if sel=='a':
#		minsel=raw_input('schiebe mehrere datensets\n waehle set min:')
#		maxsel=raw_input('waehle set max: ')
#		logtau=raw_input('um wieviele dekaden sollen die daten geschoben werden? ')
#		for i in range(int(minsel),int(maxsel)):
#			taus[i]=taus[i]+float(logtau)
#			ax.lines[i].set_xdata([brlx*10**taus[i] for brlx in brlxs[i]])
#		plt.draw()
#
#while True:
#	
#	tauout=open('tau.dat','w')
#	for i in range(0,taus.__len__()):
#		tauout.write(str(temps[i])+' '+str(taus[i])+'\n')
#	break
##		delta=[]
##		for i in range(minsel,maxsel):
##			delta.append(10**taus[i]-10**taus[minsel])
##		for i in range
#####zeichne die kurve neu
#ax.cla()
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_title('Masterkurve')
##print [brlx*10**taus[1] for brlx in brlxs[i]]
##print chis[i]
#for i in range(0,brlxs.__len__()):
#	ax.plot([brlx*10**taus[i] for brlx in brlxs[i]],
#			chis[i],
#			label=str(temps[i])+' K')
#plt.legend()
#plt.show()
#####bestimme beta
#omegataus=[]
#masterchi=[]
#for i in range(0,brlxs.__len__()):
#	for b in brlxs[i]:
#		omegataus.append(b*10**taus[i])
#for chi in chis:
#	for ch in chi:
#		masterchi.append(ch)
##fitax.plot(omegataus,masterchi,ls='None',marker='o')
#ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1./0.5,0.5))
#ax.autoscale()
#
#
#xmin=raw_input('jetzt die masterkurve nochmal fitten\n waehle xmin:')
#xmax=raw_input('waehle xmax:')
#for i in range(omegataus.__len__()-1,-1,-1):
#	if omegataus[i] < float(xmin) or omegataus[i] > float(xmax):
#		omegataus.pop(i)
#		masterchi.pop(i)
##print omegataus
##print masterchi
#params['beta'].vary=True
#params['logK_dd'].value=0
#params['logK_dd'].min=-2
#params['logK_dd'].max=2
#params['logtau'].value=0
#params['logtau'].min=-3
#params['logtau'].max=3
#omegataus=np.array(omegataus)
#out=minimize(residual,params,args=(omegataus,masterchi))
#result=omegataus+out.residual
#fit=residual(params,omegataus)
#print 'beta '+str(params['beta'].value)
#
#report_errors(params)
#####parameter updaten und ausgabe anpassen
#
#while True:
#	tauout=open('tau.dat','w')
#	for i in range(0,taus.__len__()):
#		taus[i]=taus[i]+params['logtau'].value
#		tauout.write(str(temps[i])+' '+str(taus[i])+'\n')
#	break
#omegataus=[]
#masterchi=[]
#
#ax.cla()
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_title('Masterkurve')
#
#for i in range(0,brlxs.__len__()):
#	for b in brlxs[i]:
#		omegataus.append(b*10**taus[i])
#	for chi in chis[i]:
#		masterchi.append(chi*K_dd/params['K_dd'].value)
#	ax.plot([brlx*10**taus[i] for brlx in brlxs[i]],
#			chis[i],
#			label=str(temps[i])+' K',
#			marker='o',linestyle='None')
#
#ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1./params['beta'].value,params['beta'].value),label='Modell')
##ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1.,0.5))
#
#ax.autoscale()
#plt.legend()
#plt.draw()
#		
#
#
##plt.figure(4)
##masterax=plt.axes([0.1,0.1,0.85,0.85])
##masterax.yscale=('log')
##masterax.xscale=('log')
##masterax.xlabel=(r'$\omega \tau$')
##masterax.ylabel=(r'$\frac{\chi}{K_(dd)}$')
#
#mk=raw_input('ende')
##fitpars,covmat=curve_fit(#
##		Chi_dd,	
##		brlxs[1],
##		[brlx*taus[1] for brlx in brlxs[1]],
##		p0=[1e9,1,0.5],
##		maxfev=5000
##		)
##print fitpars, covmat
#
#	#plt.plot(brlxs[1],peval(brlxs[1],plsq[0]))
#	#plt.draw()
##except ValueError: print 'what the fuck??n zum beenden'
#
#####die normierten masterdaten koennen in ein file geschrieben
#####werden.
#
#
#
#for i in range(0,temps.__len__()):
#	fout=open(str(temps[i])+' K.dat','w')
#	for ii in range(0,brlxs[i].__len__()):
#		fout.write('\n'+str(brlxs[i][ii]*10**taus[i])+' '+str(chis[i][ii]))
