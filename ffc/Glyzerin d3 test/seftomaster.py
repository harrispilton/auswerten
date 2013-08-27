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


def J_cd(omega,tau,beta):
	a=omega*beta*tau
	return (np.sin(beta*np.arctan(a))/(omega*(1+a)**(beta/2)))
def Chi_dd(omega,K_dd=1,tau=1,beta=0.5):
	return omega*3*K_dd*J_cd(omega,tau,beta)
omega=np.logspace(-3,1.5,200,10)
K_dd=1e-9
beta=0.4
tau_alpha=1
plt.figure(2)
ax=plt.axes([0.1,0.15,0.8,0.8])
ax.set_xscale('log')
ax.set_yscale('log')
plt.plot(omega,Chi_dd(omega))
plt.draw()
#def schieb(val):
#	global taus
#	slide=10.0**val
#	#if int(picker.val) == 0:
#	tau_old=taus[0]
#	taus=[tau-tau_old+slide for tau in taus]
#	#else:
#	#	tau_old=taus[int(picker.val-1)]
#	#	for i in range(int(picker.val),taus.__len__()):
#	#		taus=[tau-tau_old+slide for tau in taus]
#
#
#	for i in range(int(picker.val),brlxs.__len__()):
#		brlxs[i]=[brlx*slide for brlx in brlxs[i]]
#		ax.lines[i].set_xdata(brlxs[i])
#		brlxs[i]=[brlx/slide for brlx in brlxs[i]]
#	if int(picker.val) == 0:
#		return [tau+slide for tau in taus]
#	else:
#		#for tau in reversed(taus):#range(i,taus.__len__()):
#		#		taus[ii]=slide+taus[i-1]+tau
#		for i in xrange(int(picker.val), -1, -1):
#			taus[i]=taus[i]+slide
#		return taus
	#for i in range(0,brlxs.__len__()):
	#	brlxs[i]=[brlx*taus[i] for brlx in brlxs[i]]
#	plt.draw()
#

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
title=raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\nu$ in $MHz$")
plt.ylabel(r"$\chi$ in $s^{-2}$")
plt.xscale('log')
plt.yscale('log')
axcolor = 'lightgoldenrodyellow'


#axpicker=plt.axes([0.1,0.1,0.25,0.02],axisbg=axcolor)
#picker=Slider(axpicker,'pick set',0,sef.__len__()-0.01,valinit=0)
#axtau_c=plt.axes([0.6,0.1,0.3,0.02],axisbg=axcolor)
#stau_c=Slider(axtau_c,'log (tau_c)',-7,2,valinit=5)
#print stau_c.on_changed(schieb)##schiebe die daten durch die gegend
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
####kopplungskonstante bestimmt werden. diese muessen
####vom benutzer gewaehlt werden.




####die normierten masterdaten koennen in ein file geschrieben
####werden.
for i in range(0,temps.__len__()):
	fout=open(str(temps[i])+'.dat','w')
	for ii in range(0,brlxs[i].__len__()):
		fout.write('\n'+str(brlxs[i][ii]*10**taus[i])+' '+str(chis[i][ii]))
