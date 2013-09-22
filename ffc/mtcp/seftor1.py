#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

#text.usetex: True
###
### variablen zuweisen ordner durchfilzen
### 

sef=glob.glob('*K.sef')
sef.sort()
sdf=glob.glob('*K.sdf')
sdf.sort()
plt.figure(1)
title=raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\nu$ in $MHz$")
plt.ylabel(r"$R_1$ in $s^{-1}$")
plt.xscale('log')
plt.yscale('log')

markers=itertools.cycle(['o','s','v','x'])

sefdata=[]
temps=[]
for filename in sef:
	fin=open(filename,'r')
	sefdata=fin.readlines()
	for i in range(0,4): sefdata.pop(0)
	brlx=[]
	r1=[]
	percerr=[]
	zone=[]
	relativefile=[]
	for data in sefdata: 
		liste=data.split()
	#	liste = re.findall(r"[\w.][\f]+",data)
		brlx.append(float(liste[0]))
		r1.append(float(liste[2]))
		percerr.append(float(liste[3]))
		zone.append(int(liste[5]))
		relativefile.append(liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	print 'filename: '+filename
	#print 'zone: '+str(zone)
	print 'relativefile: '+relativefile[1]
	#print 'zone[sef.index(filename)]' + str(zone[sef.index(filename)])
	#print 'ZONE=\t'+str(zone[-2])
	temp=sdfdata[
			sdfdata.index(
				'ZONE=\t'+str(zone[
					1])+'\r\n')+7]
	temp=temp[6:]
	temp=temp.rstrip()
	temps.append(float(temp))
	#print repr(temp)
	for i,b in enumerate(brlx):
		brlx[i]=brlx[i]*1e6
	plt.plot(brlx,r1,
			label=temp+' K',
			marker=markers.next(),linestyle='None')

plt.legend()
plt.show()

