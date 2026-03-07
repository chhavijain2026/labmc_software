#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc
set dat = dataS2.dat

set dat2 = dataS2.noX.dat
awk '{ $NF=1; print }' $dat > $dat2

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0 1 2)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

   $EXEC -Di$dat2 -d1 -M1000000/$dn/200/100 -V -R$seed \
	  -Pa0/0/-5/5/0/1000/-30/30 \
          > out.S2f.$seed.dat
    sleep 1
end


# output format
# iteration#, simple_chi2, full_chi2, p1, E1 (J/mol),
# V1 (m^3/mol), A1
