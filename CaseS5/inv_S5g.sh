#!/bin/bash

#set up execution
EXEC=../src/labmc
dat=dataS5.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
for seed in 0 1 2
do
    if [ $seed -eq 0 ]
    then
       dn=1
    else
       dn=100
    fi 

    mpirun $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	   -Pf0/0/-1/3/1/5/0/0/0/0 \
	   > out.S5g.$seed.dat
    sleep 1
done

# output format
# iteration#, simple_chi2, full_chi2, p5, n5, A5
