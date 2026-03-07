#!/bin/csh -f

#set up execution - for parallelized labmc with MPI
set EXEC = ../src/labmc
set dat = path_to_data_folder/inputdata.dat

#foreach seed (0 1 2 3 4 5 6 7 8 9) # 10 parallel simulations with different initial guesses
foreach seed (0 1 2)
    if ($seed == 0) then
       set dn = 1
    else
       set dn = 100
    endif 

    mpirun -np 4 $EXEC -C3 -Di$dat -d1 -M1000000/$dn/200/100 -V -R$seed \
	  -Pa0/0/-5/5/100/1000/-10/30 \
	  -Pb0/0/-5/5/100/1000/-0/20 \
          > out.sample.$seed.dat

    sleep 1

end


# output format
# iteration#, simple_chi2, full_chi2, p1, E1, V1, n3, E3, V3, A1, A3

# dn = interval at which output is printed. To calculate ACF, \
# it is important to have all MCMC outputs printed, i.e., dn=1.\
# However, once a simulation is found to converge, subsequent \
# parallel runs, if any, may have outputs at dn>1 to reduce file size.

# The maximum number of processors = the Gibb's sample size (200 in this case).
