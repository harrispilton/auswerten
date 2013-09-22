#!/usr/bin/python
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
fig=plt.figure(1)
ax=plt.axes([0.1,0.2,0.55,0.8])
#ax2=plt.axes([0.15,0.25,0.25,0.35])
def Chi(omega,tau_c):
	return omega*tau_c/(1.+omega**2 * tau_c **2)
def update(val):
	l1=10.0**val
	print val
	print ax.lines.pop(0)
	plt.plot(omega,Chi(omega,l1))
#	for brlxs in brlx: brlxs+l1
	plt.draw()
def reset(event):
	stau_c.reset()
	
axcolor = 'lightgoldenrodyellow'
axtau_c=plt.axes([0.425,0.1,0.3,0.02],axisbg=axcolor)
stau_c=Slider(axtau_c,'log (tau_c)',-9,9,valinit=np.log10(10))
stau_c.on_changed(update)
resetax=plt.axes([0.9,0.001,0.1,0.02])
button = Button(resetax,'reset',color=axcolor,hovercolor='.975')
button.on_clicked(reset)
plt.axes([0.1,0.2,0.55,0.8])
plt.xscale('log')
plt.yscale('log')
omega=np.logspace(2,8,15,5)
plt.plot(omega,Chi(omega,2))
plt.show()
