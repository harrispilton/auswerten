#!/usr/bin/python
import glob
import re
import itertools

ordner=glob.glob('*')
ordner.sort()
ordner.pop()
print ordner
files=[]
for o in ordner:
	files.append(glob.glob(o+'/acqus')[0])
print files
#files=['3/acqus']
info=[]
for fil in files:
	index=[]
	with open(fil,'r') as fin:
		lines=fin.readlines()
		matching=[line for line in lines if r'##$NS='in line]
		index.append(lines.index(matching[0]))
		matching=[line for line in lines if r'##$TD='in line]
		index.append(lines.index(matching[0]))
		
		matching=[line for line in lines if r'##$P='in line]
		index.append(lines.index(matching[0])+1)
		matching=lines[index[-1]]#.split()##optional
		matching=[line for line in lines if r'##$D='in line]
		index.append(lines.index(matching[0])+1)
		matching=lines[index[-1]]#.split()##optional
		info.append(fil+'\n\tnumber of scans '+lines[index[0]]+'\ttime '+lines[index[1]]+'\tpulse in mus '+lines[index[2]]+'\tdelays in s '+lines[index[3]]+'\n')
print info
with open('params.dat','w') as fout:
	for inf in info:
		fout.write(inf)

