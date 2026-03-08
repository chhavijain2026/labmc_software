#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc
set dat = input/dataS7.dat

set newRuns = input/newrunsS7.dat
awk 'BEGIN { for(i=1;i<=4;i++) for(j=1;j<=8;j++) print i }' > $newRuns

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0) # 1 2)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

    $EXEC -C3 -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	  -Pa0/0/-5/5/0/0/0/0 \
	  -Pb0/0/-5/5/0/0/0/0 \
	  -B$newRuns \
          > output/out.S7fs.$seed.dat

    sleep 1

end


# output format
# iteration#, simple_chi2, full_chi2, p1, n3, A1, A3
