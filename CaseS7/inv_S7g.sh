#!/bin/csh -f

#set up execution
set EXEC = ../src/labmc
set dat = input/dataS7.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9)
foreach seed (0)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

    $EXEC -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	  -Pf0/0/-3/3/0/5/0/0/0/0 \
          > output/out.S7g.$seed.dat

    sleep 1
    
end


# output format
# iteration#, simple_chi2, full_chi2, p5, n5, A5
