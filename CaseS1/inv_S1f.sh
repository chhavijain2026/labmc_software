#!/bin/bash

#set up execution
EXEC=../src/labmc
dat=input/dataS1.dat

dat2=input/dataS1.noX.dat
awk '{$NF=1; print }' $dat > $dat2

#foreach seed (0 1 2 3 4 5 6 7 8 9)
for seed in 0 1 2
do
    if [ $seed -eq 0 ]
    then
       dn=1
    else
       dn=100
    fi 

    $EXEC -Di$dat2 -d1 -M1000000/$dn/200/100 -V -R$seed \
	   -Pa0/0/-5/5/0/0/0/0 \
	   > output/out.S1f.$seed.dat
    sleep 1
done


# output format
# iteration#, simple_chi2, full_chi2, p1, A1
