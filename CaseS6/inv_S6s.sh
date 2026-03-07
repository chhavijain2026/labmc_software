#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc 
set dataorig = dataS6.dat
set dat = dataS6.noX.dat
#set InFn = out.S6s.0.dat

awk '{$NF=1; print }' $dataorig > $dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0 1 2) #0
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

    $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	  -Pb0/0/1/5/0/1000/-30/30 \
          > out.S6s.$seed.dat

#    $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
#	  -Pb0/0/1/5/0/1000/-30/30 \
#          -I$InFn >> out.S6s.$seed.dat

    sleep 1
end


# output format
# iteration#, simple_chi2, full_chi2, n3, E3, V3, A3
