#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc
set dat = dataS8.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0) #1 2)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

    $EXEC -Di$dat -C3 -d1 -d2 -M1000000/$dn/200/100 -V -R$seed \
	  -Pa0/0/-5/5/0/1000/-30/30 \
	  -Pb0/0/1/10/0/1000/-30/30 \
          > out.S8fsX.$seed.dat
    sleep 1
end


# output format
# iteration#, simple_chi2, full_chi2, p1, E1 (J/mol), V1 (m^3/mol), 
# n3, E3, V3, X1, X2, X3, X4, X5, A1, A3
