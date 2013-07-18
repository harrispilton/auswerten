#!/usr/bin/python
import glob
import re
import numpy as np
mu_0 =1.2566e-6
h_quer = 6.626e-34/(2.*np.pi)
gamma_H=2.675e8#/2/np.pi
N_a=6.022e23
n_H=5
rho=rho_glyzerin_d3=1.26*1e6
M=M_glycerin_d3=95.0866
N=n_H*N_a*rho/M
B1=np.pi/30.*(1.+4.*(2.**0.5))*(mu_0/4./np.pi * h_quer * gamma_H **2)**2 * N
print 'B von Glycerin D3' + str(B1)
	#print N, omega,R1_0,D,B,R1_0-B/(D**1.5) *omega**0.5
#return R1_0-B/(D**1.5) *omega**0.5
mu_0 =1.2566e-6
h_quer = 6.626e-34/(2.*np.pi)
gamma_H=2.675e8#/2/np.pi
N_a=6.022e23
n_H=21
rho=rho_mTCP=1.15*1e6
M=M_mTCP=368.4
N=n_H*N_a*rho/M
B2=np.pi/30.*(1.+4.*(2.**0.5))*(mu_0/4./np.pi * h_quer * gamma_H **2)**2 * N
#print N, omega,R1_0,D,B,R1_0-B/(D**1.5) *omega**0.5
print 'B von mTCP'+str(B2)
print B2/B1
