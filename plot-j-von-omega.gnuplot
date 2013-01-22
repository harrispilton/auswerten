#!/usr/local/bin/gnuplot/
#####
##### Wir zeichnen Zunaechst die Suszeptibilitaet Chi zwei strich
##### nach der Diplomarbeit Axel Herrmann Seite 12
##### Anmerkung: leider fehlt die Normierung in der Diplomarbeit
#####
set terminal x11 persist
#set log
set title "chi for different tau_c (using a.u.)"
omega = 2.0
K_CSA(omega) = 191.3
tau_c = 103.3
delta_CSA = 123.3
delta_sigma_CSA = 12
omega_H(omega) = 18*omega
omega_P(omega) = 1.3*omega
beta=1.
h_quer = 10.3
mu_0 = 12.3
gamma_H = 11.3 
gamma_P = 100.4
r_IS = 1.1
J(omega,tau_c)=tau_c/(1. + omega **2 * tau_c **2)**beta
Chi(omega)=omega*J(omega,tau_c)
set autoscale
set xrange [1000:86200000]
set log
set xlabel "\omega _L"
set ylabel "\\omega _L * J(\omega _L)"
plot tau_c=0.5, Chi(x), tau_c=1, Chi(x), tau_c=2, Chi(x)
#pause 5
#set terminal svg 
#set output "Chi(omega).ps"
#replot
###### 
###### Nun zeichnen wir T_1 nach Adichtchev zunaechst Phosphor
######
set terminal x11 persist
set title "T_1 for different omega_L (using a.u.)"
delta_sigma_CSA(omega) = 1.#/2./omega * delta_CSA
K_CSA(omega) = 2./15. * (delta_sigma_CSA(omega) * omega)**2
T_2_CSA(omega,tau_c) = 6./(K_CSA(omega) * (3.*J(omega,tau_c) +4.* J(0,tau_c)))
T_1_CSA(omega,tau_c) = 1./(K_CSA(omega) * J(omega,tau_c))
delta_sigma_CSA(omega) = 3./2./omega * delta_CSA
#T_1_HP(omega,tau_c) = 1./(3./10. * (mu_0 / (4. * pi) * gamma_H * gamma_P * h_quer/r_IS) ** 2 * (1./3. * J(omega_H - omega_P,tau_c) + J(omega_P,tau_c) + 2*J(omega_H + omega_P,tau_c)))
T_1_HP(omega,tau_c) =  1/(1./3.*J(omega_H(omega) - omega_P(omega),tau_c) + J(omega_P(omega),tau_c) + 2. *J(omega_H(omega) + omega_P(omega),tau_c))
## Die gewuenschte Funktion
T_1(omega,tau_c) = 1./(1./T_1_HP(omega,tau_c)+1./T_1_CSA(omega,tau_c))
set autoscale
set title 'DD and CSA controlled with beta_CD =0.75 (CSA dominates)'
set log
set xrange [0.001:100]
set xlabel "log tau_c"
set ylabel "log (T_1)"
plot T_1(1,x), T_1(4,x), T_1(10,x)
set terminal svg

set output "T1 gesamt.svg" 
replot
#plot 1/T_1_CSA(1,x), 1/T_1_CSA(2,x), 1/T_1_CSA(5,x)
###
### Wir wollen 1/T_1 und betrachten den Fall, dass die
### Kopplungskonstante abhaengt von omega quadrat.(chemical shift anisotropy)
### wir betrachten hier den fall von phosphor nmr und gehen nach dem adishchev paper vor das ergebnis laesst sich leicht auf c13 uebertragen
###

## ohne werte alt:
set terminal x11 persist
set log
set title '1/T_1 mit CSA'
set ylabel '1/T_1'
set xlabel 'omega'
set xrange [.01:10]
delta_dd =1.
delta_sigma_CSA=1
K_DD=8./3.*pi**2*delta_dd**2
omega_L=8620000.
K_CSA(omega) = 2./15. * (omega * delta_sigma_CSA)**2
J(omega,tau_c)=tau_c/(1.+omega**2. * tau_c**2.)**beta
R_1(omega,tau_c) = (K_CSA(omega)*J(omega,tau_c)+K_DD*J(omega,tau_c))
# Chi_ges(omega,tau_c) = 1./(K_CSA(omega)*T_1(omega,tau_c))

print beta, delta_sigma_CSA(omega)
plot R_1(x,.001), R_1(x,4), R_1(x,10), delta_sigma_CSA=.0006 R_1(x,1)#, K_DD=1.0 R_1(x,1), R_1(x,4), K_DD=10.0 R_1(x,1), R_1(x,4)
set terminal svg
set output "Rate mit CSA.svg"
replot
pause -1

##neu mit werten
## die werte fuer omega sind von http://www-lcs.ensicaen.fr/pyPulsar/index.php/List_of_NMR_isotopes
## entnommen. die umrechnung findet sich im larmor frequenzen.ods
## alle weiteren werte sind aus dem adishchev paper
##

set terminal x11 persist
set log
set title '1/T_1 mit CSA'
set ylabel '1/T_1'
set xlabel 'omega'
set xrange [1000:8620]
delta_dd =3650.*2.*pi
delta_sigma_CSA=226.E-06 *100
K_DD=8./3.*pi**2*delta_dd**2
omega_L=8620000
K_CSA(omega) = 2./15. * (omega * delta_sigma_CSA)**2
J(omega,tau_c)=tau_c/(1.+omega**2. * tau_c**2.)**beta
R_1(omega,tau_c) = (K_CSA(omega)*J(omega,tau_c)+K_DD*J(omega,tau_c))
# Chi_ges(omega,tau_c) = 1./(K_CSA(omega)*T_1(omega,tau_c))

print beta, delta_sigma_CSA(omega)
plot R_1(x,.001), R_1(x,4), R_1(x,10), delta_sigma_CSA=.0006 R_1(x,1)#, K_DD=1.0 R_1(x,1), R_1(x,4), K_DD=10.0 R_1(x,1), R_1(x,4)
set terminal svg
set output "Rate mit CSA.svg"
replot

