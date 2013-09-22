#!/bin/bash
# sampleskript
# arrays and iteration
names=( Jennifer Tonya Anna Sadie )
for name in ${names[@]} do
echo $name
# other stuff on $name
done
for (( i = 0 ; i < ${#names[@]} ; i++ )) do
echo ${names[$i]}
# yadda yadda
done
