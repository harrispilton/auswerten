#!/bin/bash
clear
echo "delete the header of sef files. make sure there are all sef files for one measurement in the according folder"
read -p "do you want to delete all .dat files? [y] " check 
if [ "$check" == "y" ]; then 
	rm *.dat
	echo "files deleted"; 
else exit 1
fi
lsfiles="*.sef"
declare -a files=(`echo $lsfiles |sed 's/ / /g'`)
for (( i = 0; i < ${#files[@]}; i++  )) do
	cp ${files[$i]} $i.dat
	python - $i <<EOF
import sys 
import os
arg = sys.argv[1]
fin = open (sys.argv[1]+".dat","r")
lines = fin.readlines()
fin.close()
del  lines[0]
del lines[0]
del lines[0]
del lines[0]
lines.append("\n \n")
filename = lines[0]
filename = filename[-10:-6]
os.remove(sys.argv[1]+".dat")
fout = open(filename+".dat","w")
for line in lines:
        stringline=str(line)
        fout.write(stringline),
fout.close()
EOF
done
cat *.dat > all.dat
echo "${#files[@]} temperatures found"
echo "we got a all.dat file, and want to plot it"
gnuplot <<EOF
#!/usr/local/bin/gnuplot 
 set terminal x11 persist
 set palette gray
 set key left
# plot '243K.dat' using 1:3
# pause 5
tau_c=0.1
J(omega,tau_c)=tau_c/(1. + omega **2 *tau_c **2)
Chi(omega)=omega*J(omega,tau_c)
fit Chi (x) 'all.dat' index 4  using ( \$1 * 10**6):( \$2>0.002 ? \$1 * 10**6 * \$3 :1/0) via tau_c
 set log
plot Chi(x) for [i=0:${#files[@]}-1] 'all.dat' index i using ( \$1 * 10**6):( \$2>0.002 ? \$1 * 10**6 * \$3 :1/0)
#disp("blajjhg")
 pause mouse "press any key ${#files[@]}"
EOF
