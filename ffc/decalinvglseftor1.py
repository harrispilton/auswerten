#!/usr/bin/python
import glob
import re
import itertools
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

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


#text.usetex: True
###
### variablen zuweisen ordner durchfilzen
### 

sef=glob.glob('cistransdecalin/*K.sef')
sef.sort()
sdf=glob.glob('cistransdecalin/*K.sdf')
sdf.sort()
plt.ion()
plt.figure(1)
title='decalin'#raw_input("enter plot title: ")
plt.title(title)
plt.xlabel(r"$\nu$ in $MHz$")
plt.ylabel(r"$R_1$ in $s^{-1}$")
plt.xscale('log')
plt.yscale('log')
markers=get_markers()
colors=get_colors()
acolor=colors.next()

sefdata=[]
temps=[]
for filename in sef:
	amarker=markers.next()
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
		if liste[0]=='#' or float(liste[3])>40:
			pass
		else:
			brlx.append(float(liste[0]))
			r1.append(float(liste[2]))
			percerr.append(float(liste[3]))
			zone.append(int(liste[5]))
			relativefile.append('cistransdecalin/'+liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	#print 'filename: '+filename
	#print 'zone: '+str(zone)
	#print 'relativefile: '+relativefile[1]
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
	plt.scatter(brlx,r1,
			label=temp+' K cis/trans',
			marker=amarker, color='None',edgecolor=acolor)

plt.legend()
plt.show()

sef=glob.glob('transdecalin/*K.sef')
sef.sort()
sdf=glob.glob('transdecalin/*K.sdf')
sdf.sort()
markers=get_markers()
acolor=colors.next()

sefdata=[]
temps=[]
for filename in sef:
	amarker=markers.next()
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
		if liste[0]=='#' or float(liste[3])>40:
			pass
		else:
			brlx.append(float(liste[0]))
			r1.append(float(liste[2]))
			percerr.append(float(liste[3]))
			zone.append(int(liste[5]))
			relativefile.append('transdecalin/'+liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	#print 'filename: '+filename
	#print 'zone: '+str(zone)
	#print 'relativefile: '+relativefile[1]
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
	plt.scatter(brlx,r1,
			label=temp+' K trans',
			marker=amarker, color='None',edgecolor=acolor)

plt.legend()
plt.show()
sef=glob.glob('decalin in otp/*K.sef')
sef.sort()
sdf=glob.glob('decalin in otp/*K.sdf')
sdf.sort()
markers=get_markers()
acolor=colors.next()

sefdata=[]
temps=[]
for filename in sef:
	amarker=markers.next()
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
		if liste[0]=='#' or float(liste[3])>40:
			pass
		else:
			brlx.append(float(liste[0]))
			r1.append(float(liste[2]))
			percerr.append(float(liste[3]))
			zone.append(int(liste[5]))
			relativefile.append('decalin in otp/'+liste[6])
	fin2=open(relativefile[1],'r')
	sdfdata=fin2.readlines()
	#print 'filename: '+filename
	#print 'zone: '+str(zone)
	#print 'relativefile: '+relativefile[1]
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
	plt.scatter(brlx,r1,
			label=temp+' K in otp',
			marker=amarker, color='None',edgecolor=acolor)

plt.legend()
plt.show()
sdflk=raw_input('ende')
