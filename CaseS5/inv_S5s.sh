#!/bin/bash

#set up execution
EXEC=../src/labmc
dat=input/dataS5.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
for seed in 0 1 2
do
    if [ $seed -eq 0 ]
    then
       dn=1
    else
       dn=100
    fi 

    $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	   -Pb0/0/1/5/0/0/0/0 \
	   > output/out.S5s.$seed.dat
    sleep 1
done

wait

# output format
# iteration#, simple_chi2, full_chi2, n3, A3
