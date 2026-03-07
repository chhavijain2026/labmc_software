#!/bin/bash

#set up execution
EXEC=../src/labmc
datorig=dataS4.dat
dat=dataS4.noX.dat

awk '{$NF=1; print }' $datorig > $dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
for seed in 0 #1 2
do
    if [ $seed -eq 0 ]
    then
       dn=1
    else
       dn=100
    fi 

    $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	   -Pa0/0/-5/5/0/1000/0/0 \
	   > out.S4f.$seed.dat
    sleep 1
done


# output format
# iteration#, simple_chi2, full_chi2,  p1, E1 (J/mol), A1
