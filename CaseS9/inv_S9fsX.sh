#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc
set dat = input/dataS9.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0 1)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

   mpirun $EXEC -Di$dat -C5 -d1 -d2 -M1000000/$dn/300/100 -V -R$seed \
	  -Pd0/0/1/10/-2/2/0/1000/-30/30 \
	  -Pc0/0/-5/5/-2/2/0/1000/-30/30 \
          > output/out.S9fsX.$seed.dat
    sleep 1
end


# output format
# iteration#, simple_chi2, full_chi2, p2, r2, E2 (J/mol), V2 (m^3/mol), 
# r4, n4, E4, V4, X1, X2, ..., X12, A2, A4
