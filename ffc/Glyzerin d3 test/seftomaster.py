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
from lmfit import  minimize, Parameters, report_errors

def residual(params,xdata,ydata=None,releps=None):
	K_dd=params['K_dd'].value
	tau=params['tau'].value
	beta=params['beta'].value
	if ydata is None:
		return Chi_dd
	if releps is None:
		return (ydata-Chi_dd(xdata,K_dd,tau,beta))
	return ((ydata-Chi_dd(xdata,K_dd,tau,beta))*releps)
def J_p(omega,tau):
	return (tau/(1+(omega*tau)**2))
def J_cd(omega,tau,beta):
	a=omega*beta*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1+a)**(beta/2)))
def Chi_dd(omega,K_dd=1e-8,tau=1e-6,beta=0.5):
	return omega*3*K_dd*J_cd(omega,tau,beta)
omega=np.logspace(-3,1.5,200,10)
#K_dd=1e-9
#beta=0.4
#tau_alpha=1
params= Parameters()
params.add('logK_dd',value=8.0,min=7,max=10)
params.add('K_dd',expr='(10.0**logK_dd)')
params.add('logtau',value=-6.0,min=-10,max=10)
params.add('tau',expr='(10.0**logtau)')
params.add('beta',value=0.45,vary=False)


plt.figure(2)
ax=plt.axes([0.1,0.15,0.8,0.8])
ax.set_xscale('log')
ax.set_yscale('log')
plt.plot(omega,Chi_dd(omega))
plt.draw()
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
title='dummer fit'#raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\nu$ in $MHz$")
plt.ylabel(r"$\chi$ in $s^{-2}$")
plt.xscale('log')
plt.yscale('log')
axcolor = 'lightgoldenrodyellow'

markers=itertools.cycle(['o','s','v','x'])

sefdata=[]
temps=[]
brlxs=[]
chis=[]
omegas=[]
taus=[]##liste mit log10(tau_strukturrelaxation)
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	brlx=[]
	chi=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		brlx.append(float(liste[0])*1e6)
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
	chis.append(chi)
	taus.append(0.0)
	#print repr(temp)
	ax.plot(brlx,chi,
			label=temp+' K',
			marker=markers.next(),linestyle='None')

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
ax.legend()
plt.show()
#while True:
#	selecttofit=raw_input("waehle die datensets mit alphapeak (0-"+str(temps.__len__())+") abbrechen mit n:  ")
#	if selecttofit==n:break
#	brlxs[int(selecttofit)]
####im folgenden kann man einzelne daten per eingabedialog verschieben
####die strukturrelaxationszeiten werden logarithmisiert in
####einer liste abgelegt
while True:
	sel=raw_input("Waehle den datenset: ")
	try: 
		int(sel)
		print "aktuelles log tau: "+str(taus[int(sel)])
		while True:
			logtau=raw_input("neues tau: ")
			if logtau=='n':break
			tau=10**float(logtau)
			taus[int(sel)]=float(logtau)
			ax.lines[int(sel)].set_xdata([brlx*tau for brlx in brlxs[int(sel)]])
			plt.draw()
	except ValueError: print 'n zum beenden'
	if sel=='n':break


####aus datensaetzen mit peak im frequenzfenster kann die
####kopplungskonstante bestimmt werden. datensaetze muessen
####vom benutzer gewaehlt werden.

#k=raw_input("Schaetze die Kopplungskonstante:  ") 
#try:
	#k=float(k)
seti=[]
while True:
	seti.append(raw_input("waehle ein datenset mit maximum: "))
	try: int(seti[-1])
	except ValueError: print 'n zum beenden'
for i in seti:
	brlxs[i]=np.array(brlxs[i])
	k,tau,beta=1e-5,2e-6,0.5
	out = minimize(residual, params,args=(brlxs[i],chis[i]),method=('leastsq'))
	result=brlxs[i]+out.residual
	fit = residual(params,brlxs[i])
	print fit
	print result
	print params['K_dd'].value
	print params['beta'].value
	print params['tau'].value
	
	report_errors(params)
	plt.figure(3)
	fitax=plt.axes([0.1,0.1,0.8,0.8])
	fitax.set_xscale('log')
	fitax.set_yscale('log')
	plt.plot(brlxs[i],Chi_dd(brlxs[i],params['K_dd'].value,params['tau'].value,0.5))
	plt.plot(brlxs[i],chis[i])
	plt.show()
	k=raw_input('fuer naechster fit ente druecken')
k=raw_input('ende')
#fitpars,covmat=curve_fit(#
#		Chi_dd,	
#		brlxs[1],
#		[brlx*taus[1] for brlx in brlxs[1]],
#		p0=[1e9,1,0.5],
#		maxfev=5000
#		)
#print fitpars, covmat

	#plt.plot(brlxs[1],peval(brlxs[1],plsq[0]))
	#plt.draw()
#except ValueError: print 'what the fuck??n zum beenden'

####die normierten masterdaten koennen in ein file geschrieben
####werden.


for i in range(0,temps.__len__()):
	fout=open(str(temps[i])+' K.dat','w')
	for ii in range(0,brlxs[i].__len__()):
		fout.write('\n'+str(brlxs[i][ii]*10**taus[i])+' '+str(chis[i][ii]))
