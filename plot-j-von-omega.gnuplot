#!/usr/local/bin/gnuplot/
#####
##### Wir zeichnen Zunaechst die Suszeptibilitaet Chi zwei strich
##### nach der Diplomarbeit Axel Herrmann Seite 12
##### Anmerkung: leider fehlt die Normierung in der Diplomarbeit
#####
set terminal x11 persist
#set log
set title "chi for different tau_c (using a.u.)"
K_CSA(omega) = 191.3
tau_c = 103.3
delta_CSA = 123.3
delta_sigma_CSA = 12
omega_H = 1.3
omega_P = 1.3
beta=1.
h_quer = 10.3
mu_0 = 12.3
gamma_H = 11.3 
gamma_P = 100.4
r_IS = 1.1
J(omega,tau_c)=tau_c/(1. + omega **2 * tau_c **2)**beta
Chi(omega)=omega*J(omega,tau_c)
set autoscale
set xrange [0.01:100]
set log
set xlabel "\omega _L"
set ylabel "\\omega _L * J(\omega _L)"
plot tau_c=0.5, Chi(x), tau_c=1, Chi(x), tau_c=2, Chi(x)
#set terminal svg 
#set output "Chi(omega).ps"
#replot
###### 
###### Nun zeichnen wir T_1 nach Adichtchev zunaechst Phosphor
######
set terminal x11 persist
set title "T_1 for different omega_L (using a.u.)"
delta_sigma_CSA(omega) = 1.#/2./omega * delta_CSA
K_CSA(omega) =10*omega**2#2./15. * (delta_sigma_CSA(omega) * omega)**2
T_2_CSA(omega,tau_c) = 6./(K_CSA(omega) * (3.*J(omega,tau_c) +4.* J(0,tau_c)))
T_1_CSA(omega,tau_c) = 1./(K_CSA(omega) * J(omega,tau_c))
delta_sigma_CSA(omega) = 3./2./omega * delta_CSA
#T_1_HP(omega,tau_c) = 1./(3./10. * (mu_0 / (4. * pi) * gamma_H * gamma_P * h_quer/r_IS) ** 2 * (1./3. * J(omega_H - omega_P,tau_c) + J(omega_P,tau_c) + 2*J(omega_H + omega_P,tau_c)))
T_1_HP(tau_c) =  1/(J(0.1,tau_c) + J(2,tau_c) + J(4,tau_c))
## Die gewuenschte Funktion
T_1(omega,tau_c) = 1./(1./T_1_HP(tau_c))
set autoscale
set title 'DD controlled'
set log
set xrange [0.001:100]
set xlabel "log tau_c"
set ylabel "log (T_1)"
plot T_1(1,x), T_1(10,x), T_1(100,x)
set terminal svg

set output "T1 HP.ps"
replot
#plot 1/T_1_CSA(1,x), 1/T_1_CSA(2,x), 1/T_1_CSA(5,x)
