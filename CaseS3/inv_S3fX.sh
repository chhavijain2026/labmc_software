#!/bin/bash

#set up execution
EXEC=../src/labmc
dat=input/dataS3.dat

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
	   -Pa0/0/-5/5/0/0/0/0 \
	   > output/out.S3fX.$seed.dat
    sleep 1
done


# output format
# iteration#, simple_chi2, full_chi2, p1, X1, X2, X3, A1
