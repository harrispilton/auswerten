#!/usr/local/bin/gnuplot/
##hopefully quick and not too dirty
set terminal x11 persist
set log
set title 'NMR suszeptibilitaet'
set xlabel 'omega'
set ylabel 'chi 2strich'
set key outside
plot 	'230K.dat' using ($1 * 10**6):( $2>0.002 ? $1 * 10**6 *  $3 :1/0) title '230', \
	'330K.dat' using ($1 * 10**6):( $2>0.002 ? $1 * 10**6 *  $3 :1/0) title '330'

set terminal svg
set output 'chibrut.svg'
replot
