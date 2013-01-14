#!/usr/local/bin/gnuplot/
##hopefully quick and not too dirty
set terminal x11 persist
set log
set title 'NMR suszeptibilitaet'
set xlabel 'omega'
set ylabel 'chi 2strich'

plot 'all.dat' using ($1 * 10**6):( $2>0.002 ? $1 * 10**6 *  $3 :1/0)
set terminal svg
set output 'chibrut.svg'
replot
