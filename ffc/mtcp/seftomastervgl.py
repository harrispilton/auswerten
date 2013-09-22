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
#,,from lmfit import  minimize, Parameters, report_errors

def residual(params,xdata,ydata=None,releps=None):
	K_dd=params['K_dd'].value
	tau=params['tau'].value
	beta=params['beta'].value
	if ydata is None:
		return Chi_dd
	if releps is None:
		return (ydata-Chi_dd(xdata,K_dd,tau,beta))
	return ((ydata-Chi_dd(xdata,K_dd,tau,beta))*releps)
def J_p(nu,tau):
	omega=nu*2*np.pi
	return (tau/(1+(omega*tau)**2))
def J_cd(nu,tau,beta):
	omega=nu*2*np.pi
	a=omega*beta*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1+a)**(beta/2)))
def Chi_dd(nu,K_dd=1,tau=1,beta=0.5):
	omega=nu*2*np.pi
	return omega*3*K_dd*J_cd(nu,tau,beta)

#K_dd=1e-9
#beta=0.4
#tau_alpha=1
#omega=np.logspace(-3,1.5,200,10)
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
c=[]
for i in np.arange(40):c.append(cm.jet(i/40.))
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

markers=itertools.cycle(['o','s','v','x'])

sefdata=[]
temps=[]
brlxs=[]
chis=[]
omegas=[]
percerrs=[]
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
	percerrs.append(percerr)
	#print repr(temp)
	ax.plot(brlx,chi,
			label=temp+' K',
			marker=markers.next(),linestyle='None')

for i in range(0, temps.__len__()):print str(i)+':   ', str(temps[i])
ax.legend()
plt.draw()
#while True:
#	selecttofit=raw_input("waehle die datensets mit alphapeak (0-"+str(temps.__len__())+") abbrechen mit n:  ")
#	if selecttofit==n:break
#	brlxs[int(selecttofit)]
####Normierung der Hoehe:
####aus datensaetzen mit peak im frequenzfenster kann die
####kopplungskonstante bestimmt werden. datensaetze muessen
####vom benutzer gewaehlt werden.

#k=raw_input("Schaetze die Kopplungskonstante:  ") 
#try:
	#k=float(k)
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
#	out = minimize(residual, params,args=(brlxs[i],chis[i],percerrs[i]),method=('leastsq'))
#	result=brlxs[i]+out.residual
#	fit = residual(params,brlxs[i])
#	print fit
#	print result
#	ks.append(params['K_dd'].value)
#	print params['beta'].value
#	print params['tau'].value
#	
#	report_errors(params)
#	plt.figure(3)
#
#	fitax=plt.axes([0.1,0.1,0.8,0.8])
#	fitax.set_xscale('log')
#	fitax.set_yscale('log')
#	plt.plot(brlxs[i],Chi_dd(brlxs[i],params['K_dd'].value,params['tau'].value,0.5))
#	plt.plot(brlxs[i],chis[i])
#	plt.show()

#	mk=raw_input('fuer naechster fit ente druecken')
#K_dd=0.0
K_dd=5.52e8
#for k in ks:
#	K_dd=K_dd+k
#K_dd=K_dd/ks.__len__()
#for i in range(0,chis.__len__()):
#	chis[i]=[chi/K_dd for chi in chis[i]]
#	ax.lines[i].set_ydata(chis[i])	
#om=np.logspace(-6,3,80)
#ax.plot(om,Chi_dd(om,1.,1.),label='fit')
#plt.draw()
#while True:
#	print 'K alt='+str(K_dd)
#	kneu=raw_input('neues K: ')
#	if kneu=='n': break
#	for i in range(0,chis.__len__()):
#		chis[i]=[chi/float(kneu)*K_dd for chi in chis[i]]
#		ax.lines[i].set_ydata(chis[i])	
#		plt.draw()

#####importiere die dielektrischen daten
plt.xscale('log')
plt.yscale('log')


while True:
	masterin=open('externe_daten/mtcp_master_dsX\'\'.dat','r')
	dsdat=masterin.readlines()
	omegatau=[]
	chi=[]
	for dat in dsdat:
		liste=dat.split()
		omegatau.append(float(liste[0]))
		chi.append(float(liste[1]))
	ax.plot(omegatau,chi,label='ds',marker='o')
	plt.legend()
	plt.autoscale()
	plt.draw()
	plt.show()
	break

####importiere die pcs daten
while True:
	masterin=open('externe_daten/pcs fit.dat','r')
	dsdat=masterin.readlines()
	omegatau=[]
	chi=[]
	for dat in dsdat:
		liste=dat.split()
		omegatau.append(float(liste[0]))
		chi.append(float(liste[1]))
	ax.plot(omegatau,chi,label='pcs',marker='o')
	plt.legend()
	plt.autoscale()
	plt.draw()
	plt.show()
	break

####Uebereinanderschieben der Daten:
####die taus koennen von Hand eingegben werden,
####Fits bringen hier niemanden weiter
####Im Folgenden kann man einzelne Daten per Eingabedialog verschieben
####die Strukturrelaxationszeiten werden logarithmisiert in
####einer Liste abgelegt
while True:
	tauin=open('tau.dat','r')
	lines=tauin.readlines()
	if lines.__len__()==taus.__len__():
		for i in range(0,lines.__len__()):
			liste=lines[i].split()
			taus[i]=float(liste[1])
			ax.lines[i].set_xdata([brlx*10**taus[i] for brlx in brlxs[i]])
		plt.autoscale()
		plt.draw()
	else: a=raw_input('laenge der tau stimmt nicht...')
	break

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
	if sel=='n':
		break
	if sel=='a':
		minsel=raw_input('schiebe mehrere datensets\n waehle set min:')
		maxsel=raw_input('waehle set max: ')
		logtau=raw_input('um wieviele dekaden sollen die daten geschoben werden? ')
		for i in range(int(minsel),int(maxsel)):
			taus[i]=taus[i]+float(logtau)
			ax.lines[i].set_xdata([brlx*10**taus[i] for brlx in brlxs[i]])
		ax.autoscale()
		plt.draw()
with open('tau.dat','w') as tauout:
	for i  in range(0,taus.__len__()):
		tauout.write(str(temps[i])+' '+str(taus[i])+'\n')
while True:
	
	tauout=open('tau.dat','w')
	for i in range(0,taus.__len__()):
		tauout.write(str(temps[i])+' '+str(taus[i])+'\n')
	break
#		delta=[]
#		for i in range(minsel,maxsel):
#			delta.append(10**taus[i]-10**taus[minsel])
#		for i in range
####zeichne die kurve neu
for i in range(0,taus.__len__()):
	ax.lines[i].set_ydata([chi*7.075e-10 for chi in chis[i]])
	plt.draw()
plt.show()

fout=open('master.dat','w')
for i in range(0,taus.__len__()):
	for ii in range(0,chis[i].__len__()):
		fout.write(str(brlxs[i][ii]*10**taus[i])+' '+str(chis[i][ii]*7.075e-10)+'\n')
i=raw_input('ende')
ax.cla()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Masterkurve')
#print [brlx*10**taus[1] for brlx in brlxs[i]]
#print chis[i]
for i in range(0,brlxs.__len__()):
	ax.plot([brlx*10**taus[i] for brlx in brlxs[i]],
			chis[i],
			label=str(temps[i])+' K')
	plt.draw()
plt.legend()
plt.show()
####bestimme beta
omegataus=[]
masterchi=[]
mastererr=[]
for i in range(0,brlxs.__len__()):
	for b in brlxs[i]:
		omegataus.append(b*10**taus[i])
for chi in chis:
	for ch in chi:
		masterchi.append(ch)
for percerr in percerrs:
	for p in percerr:
		mastererr.append(p)
ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1.,0.5))
ax.autoscale()
plt.draw()


xmin=raw_input('jetzt die masterkurve nochmal fitten\n waehle xmin:')
xmax=raw_input('waehle xmax:')
for i in range(omegataus.__len__()-1,-1,-1):
	if omegataus[i] < float(xmin) or omegataus[i] > float(xmax):
		omegataus.pop(i)
		masterchi.pop(i)
		mastererr.pop(i)
params['beta'].vary=True
params['logK_dd'].value=0
params['logK_dd'].min=-2
params['logK_dd'].max=2
params['logtau'].value=0
params['logtau'].min=-3
params['logtau'].max=3
omegataus=np.array(omegataus)
out=minimize(residual,params,args=(omegataus,masterchi,mastererr))
result=omegataus+out.residual
fit=residual(params,omegataus)
print 'beta '+str(params['beta'].value)
report_errors(params)

####parameter updaten und ausgabe anpassen
while True:
	tauout=open('tau.dat','w')
	for i in range(0,taus.__len__()):
		taus[i]=taus[i]+params['logtau'].value
		tauout.write(str(temps[i])+' '+str(taus[i])+'\n')
	break
omegataus=[]
masterchi=[]

ax.cla()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Masterkurve')

for i in range(0,brlxs.__len__()):
	for b in brlxs[i]:
		omegataus.append(b*10**taus[i])
	for chi in chis[i]:
		masterchi.append(chi*K_dd/params['K_dd'].value)
	ax.plot([brlx*10**taus[i] for brlx in brlxs[i]],
			chis[i],
			label=str(temps[i])+' K',
			marker='o',linestyle='None',
			color=c[int((float(temps[i])-160)/6.)]
			)
	plt.draw()
ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1.,params['beta'].value),label='Modell')
plt.draw()
try:
	while True:
		tauc=raw_input('passe tau an (faktor): ')
		if tauc != 'n':
			taus=[np.log10(float(tauc))+tau for tau in taus]
			for i in range(0,ax.lines.__len__()-1):
				ax.lines[i].set_xdata([b*10**taus[i] for b in brlxs[i]])
			xmin=float(xmin)*float(tauc)
			xmax=float(xmax)*float(tauc)
		k=raw_input('passe k an (faktor): ')
		if k!= 'n':
			params['K_dd'].value=params['K_dd'].value*float(k)
			for i  in range(0,ax.lines.__len__()-1):
				ax.lines[i].set_ydata([ch*params['K_dd'].value for ch in chis[i]])
		plt.draw()
		if tau=='n' and k == 'n': break
except ValueError: 
	print 'Value Error'


#ax.plot(sorted(omegataus),Chi_dd(np.array(sorted(omegataus)),1.,1.,0.5))

ax.xlabel=(r'$\omega \tau$')
ax.ylabel=(r'$\frac{\chi}{K_(dd)}$')
ax.autoscale()
plt.legend()
plt.draw()

####2nd scaling
kdds=[]
cbrlxs=[]
cchis=[]
cpercerrs=[]
params['tau'].vary=False
params['beta'].vary=False
print 'xmax'
print xmax
print 'xmin'
print xmin

####waehle zunaechst die wertepaare die fuers fitten in frage
####kommen aus (massgeblich sind die vorher festgelegten xmin xmax)
for i in range(0,brlxs.__len__()):
	cbrlxs.append([])
	cchis.append([])
	cpercerrs.append([])
	for ii in range(0,brlxs[i].__len__()):
		if xmin<brlxs[i][ii]<xmax:
			cbrlxs[i].append(brlxs[i][ii])
			cchis[i].append(chis[i][ii])
			cpercerrs[i].append(percerrs[i][ii])
for i in range(cbrlxs.__len__()-1,-1,-1):
	if cbrlxs[i].__len__()<3:
		print temps[i]
		cbrlxs.pop(i)
		cchis.pop(i)
		cpercerrs.pop(i)
####minimiere den abstand zum model ohne die ywerte zu schieben
for i in range(0,cbrlxs.__len__()):
	cbrlxs[i]=np.array(cbrlxs[i])
	out=minimize(residual,params,args=(cbrlxs[i],cchis[i],cpercerrs[i]),method=('leastsq'))
	print report_errors(params)
	kdds.append(params['K_dd'].value)
	chis[i]=[kdds[i]*ch for ch in chis[i]]
	ax.lines[i].set_ydata(chis[i])
plt.draw()

	

#plt.figure(4)
#masterax=plt.axes([0.1,0.1,0.85,0.85])
#masterax.yscale=('log')
#masterax.xscale=('log')
mk=raw_input('ende')
with open('tau.dat','w') as fout:
	for i in range(0,taus.__len__()):
		fout.write(str(temps[i])+' '+str(taus[i])+'\n')

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



#for i in range(0,temps.__len__()):
#	fout=open(str(temps[i])+' K.dat','w')
#	for ii in range(0,brlxs[i].__len__()):
#		fout.write('\n'+str(brlxs[i][ii]*10**taus[i])+' '+str(chis[i][ii]))
