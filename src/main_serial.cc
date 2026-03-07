/*
 * Latest Revision in Jan'18
 * This code fits gradient b/w each data-pair. All pairs within a run are considered
 * Simple chi^2 is used as cost function (MKK15)
 * parallelized
 * slightly different input file options than the MKK15 code
 * new flowlaws added: 3 eqns for Peierls & GBS (all dry)
 * More new flow laws added: a 4th eqn for Peierls & wet GBS
 * Added a new option to regroup data pairs when inter-run bias=0 within a study
 *               ( to save time. CG search for ALL possible pairs takes v long)
 * Added a new option to group sticky parameters together
 * Modified sticky parameter option to allow multiple calls to -G
 *
 * usage: labmc -Di<data> -P<flowlaw_para> -S<flowlaw_seq>
 *              -M<niter>/<dn>/<m>/<dr> [ -F<i1/i2> -B<min>/<max> ]
 *              [ -R<seed> -C<n>/<dcg> -V ] [-I<inFile>] [-i]
 *              [-Do<outData file>] 
 *
 *        -D - sets data file
 *           -Di - sets input data file that has to be read
 *                 First 8 items should be (T dT p dp e de sig dsig)
 *                 Last item should be run id (must be integer)
 *           -Do - sets output file in which randomized data is printed
 *                 The format is same as the input data file except that
 *                 the first column is the data number. When this sequence
 *                 is reset, it marks the new randomized data set
 *
 *        -d1  grain size (d dd)
 *        -d2  water content (COH dCOH)
 *        -d3  oxygen fugacity (fO2 dfO2)
 *        -d4  melt fraction (phi dphi)
 *
 *        -M - sets MCMC configuration
 *             <ninter> = maximum number of MC iterations
 *             <dn>     = output interval
 *             <m>      = number of trial calculations in rejection loop
 *             <dr>     = data randomization interval
 *
 *        -P - sets a parallel flow law
 *
 *             [dry diffusion creep]
 *             -Pa<facA_min>/<facA_max>/<mmin>/<mmax>
 *                /<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [dry dislocation creep]
 *             -Pb<facA_min>/<facA_max>/<nmin>/<nmax>
 *                /<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [wet diffusion creep]
 *             -Pc<facA_min>/<facA_max>/<mmin>/<mmax>
 *                 /<rmin>/<rmax>/<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [wet dislocation creep]
 *             -Pd<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<rmin>/<rmax>/<Emin>/<Emax>/<Vmin>/<Vmax>
 *
 *             [gen flow creep]
 *             -Pe<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<smin>/<smax>/<Qmin>/<Qmax>
 *
 *             [Peierls mechanism (dry)]
 *             -Pp<facA_min>/<facA_max>/<sigP0_min>/<sigP0_max>
 *              /<E_min>/<E_max>/<V_min>/<V_max>-q<q>s<s>
 *
 *             [Peierls mechanism2 (dry)]-no sigma^2 term
 *             -Pq<facA_min>/<facA_max>/<sigP0_min>/<sigP0_max>
 *              /<E_min>/<E_max>/<V_min>/<V_max>-q<q>s<s>
 *
 *             [Peierls mechanism3 (dry)]-no pressure effect but with sigma^2
 *             -Pq<facA_min>/<facA_max>/<sigP0_min>/<sigP0_max>
 *              /<E_min>/<E_max>-q<q>s<s>
 *
 *             [Peirls mechanism4 (dry)]-no pressure effect & without sigma^2
 *             -Ps<facA_min>/<facA_max>/<sigP0_min>/<sigP0_max>
 *               /<E_min>/<E_max>-q<q>s<s>
 *
 *             [dry GBS]
 *             -Pf<facA_min>/<facA_max>/<m_min>/<m_max>
 *              /<n_min>/<n_max>/<E_min>/<E_max>/<V_min>/<V_max>
 *
 *             [wet GBS]
 *             -Pg<facA_min>/<facA_max>/<m_min>/<m_max>
 *              /<n_min>/<n_max>/<r_min>/<r_max>/<E_min>/<E_max>/<V_min>/<V_max>
 *
 *        -S - sets a sequential flow law
 *
 *             [dislocation creep]
 *             -Sd<facA_min>/<facA_max>/<nmin><nmax>
 *                 /<mmin>/<mmax>/<Emin>/<Emax>/<Vmin>/<Vmax>
 *             [gen flow creep]
 *             -Se<facA_min>/<facA_max>/<nmin>/<nmax>
 *                 /<smin>/<smax>/<Qmin>/<Qmax>
 *
 *        -F - fix param(i1) to param(i2)
 *             (i.e., two parameters always share the same value
 *              e.g., -F3/2 & -F4/2 --- param(3) and param(4) will be
 *                                      fixed to param(2)
 *              note: -F3/2 & -F4/3 probably won't give the correct result.)
 *
 *        -G - group sticky parameters (x1/../xn) together so they are 
 *             randomized together. They are randomized together after 
 *             every interval di
 *             -G<x1/x1/../xn> -- e.g., -G2/5: param(2) and 
 *                   param(5) are sampled together.
 *
 *        -B - <name of file> which contains new grouping of runs: sets new_pairs=true
 *        -R - sets seed for rand()
 *        -C - forces to use conjugate gradient search for `best-fit'
 *             scaling factors
 *             <n> = the number of different initial values to be tested
 *        -I - <name of file> whose last line will be read & will serve as 
 *                            the initial MCMC model
 *            
 *        -V - sets verbose mode
 *
 * Jun Korenaga
 * Summer 2008
 * Revised: Fall 2012 (for lambda-approach)
 * Revised: Fall 2014 (for data-randomization approach)
 * 
 * CJ - Revised in 2015 for new gradient-fitting approach, parallelized the code,
 *       and added new mechanisms
 * CJ - Revised in 2017 to incorporate sticky parameters. Also, included an option
 *      for reduced pairing in the conjugate gradient part to save time.
 * CJ - Revised in February 2018 to include provision for fixed A (statmodel modified)
 *      and modified nbias factor so that no. of inter-run bias = 0 when for only 1 run  
 */

#include <iostream>
#include <fstream>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <time.h>
#include "labmc_nr.h"
#include "constants.h"
#include "array.h"
#include "util.h"
#include "parameter.h"
#include "state.h"
#include "flowlaw.h"
#include "statmodel.h"
#include <stdio.h>
#include <cstring>
#include <chrono>

//#define _CJ_DEBUG_
/*  debug  */
//#ifdef _CJ_DEBUG_
//#include "debug.h"
//#endif /* _CJ_DEBUG_ */
#define _TIMECALC_

// ---- Serial timing helper (replacement for MPI_Wtime) ----
static double wtime()
{
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    const auto now = clock::now();
    return std::chrono::duration<double>(now - t0).count();
}

int main(int argc, char **argv)
{
    int nerror=0, start=0;
    char *dfn, *infn, *outD, *rfn;
    bool getD=false, getM=false, useCG=false, readfile=false, verbose=false;
    bool outRanData=false, new_pairs=false, sticky=false, isSticky=false;
    int max_iter, dn_out, max_m, dn_ran, nNewRuns, gn=0;
    Array1d<int> ifix, isrc, sgrp;
    int ran_seed=1;
    StatModel model;

    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'D':
	      switch(argv[i][2]){
	      case 'i':
		{
		  dfn = &argv[i][3];
		  getD = true;
		  break;
		}
	      case 'o':
		{
		  outRanData = true;
		  outD = &argv[i][3];
		  break;
		}
	      }
	    case 'd':
	    {
		int istate = atoi(&argv[i][2]);
		switch(istate){
		case 1:
		    model.stateToRead(State::grain_size);
		    break;
		case 2:
		    model.stateToRead(State::water_content);
		    break;
		case 3:
		    model.stateToRead(State::oxygen_fugacity);
		    break;
		case 4:
		    model.stateToRead(State::melt_fraction);
		    break;
		}
		break;
	    }
	    case 'M':
		if (sscanf(&argv[i][2],
			   "%d/%d/%d/%d", &max_iter, &dn_out, 
			   &max_m, &dn_ran) != 4){
		    cerr << "invalid -M option\n";
		    nerror++;
		}
		getM = true;
		break;
	    case 'P':
		switch(argv[i][2]){
		case 'a':
		{
		  double a1,a2,m1,m2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&e1,&e2,&v1,&v2) != 8){
			cerr << "invalid -Pa option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDiffDry(a1,a2,m1,m2,e1,e2,v1,v2));
		    break;
		}

		case 'b':
		{
		  double a1,a2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&e1,&e2,&v1,&v2) != 8){
			cerr << "invalid -Pb option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDisDry(a1,a2,n1,n2,e1,e2,v1,v2));
		    break;
		}

		case 'p':
		{
		  double q1,q2,a1,a2,sp1,sp2,e1,e2,v1,v2;

		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf-q%lf/s%lf",
			     &a1,&a2,&sp1,&sp2,&e1,&e2,&v1,&v2,&q1,&q2) != 10){
			cerr << "invalid -Pp option\n";
			nerror++;
		    }
		    model.addParallel(new
			     FlowLawPeierls(a1,a2,sp1,sp2,e1,e2,v1,v2,q1,q2));
		    break;
		}
		case 'q':
		{
		  double q1,q2,a1,a2,sp1,sp2,e1,e2,v1,v2;

		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf-q%lf/s%lf",
			     &a1,&a2,&sp1,&sp2,&e1,&e2,&v1,&v2,&q1,&q2) != 10){
			cerr << "invalid -Pq option\n";
			nerror++;
		    }
		    model.addParallel(new
			     FlowLawPeierls2(a1,a2,sp1,sp2,e1,e2,v1,v2,q1,q2));
		    break;
		}
		case 'r':
		{
		  double q1,q2,a1,a2,sp1,sp2,e1,e2;

		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf-q%lf/s%lf",
			     &a1,&a2,&sp1,&sp2,&e1,&e2,&q1,&q2) != 8){
			cerr << "invalid -Pr option\n";
			nerror++;
		    }
		    model.addParallel(new
			     FlowLawPeierls3(a1,a2,sp1,sp2,e1,e2,q1,q2));
		    break;
		}
		case 's':
                {
                  double q1,q2,a1,a2,sp1,sp2,e1,e2;

                    if (sscanf(&argv[i][3],
                               "%lf/%lf/%lf/%lf/%lf/%lf-q%lf/s%lf",
                             &a1,&a2,&sp1,&sp2,&e1,&e2,&q1,&q2) != 8){
                        cerr << "invalid -Ps option\n";
                        nerror++;
                    }
                    model.addParallel(new
                             FlowLawPeierls4(a1,a2,sp1,sp2,e1,e2,q1,q2));
                    break;
                }
		case 'c':
		{
		  double a1,a2,m1,m2,r1,r2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Pc option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDiffWet(a1,a2,m1,m2,r1,r2,
						     e1,e2,v1,v2));
		    break;
		}
		case 'd':
		{
		  double a1,a2,r1,r2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Pd option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawDisWet(a1,a2,n1,n2,r1,r2,
						     e1,e2,v1,v2));
		    break;
		}
		case 'e':
		{
		  double a1,a2,s1,s2,n1,n2,e1,e2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&s1,&s2,&e1,&e2) != 8){
			cerr << "invalid -Pe option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawGen(a1,a2,n1,n2,s1,s2,
						     e1,e2));
		    break;
		}
		case 'f':
		{
		  double a1,a2,m1,m2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&m1,&m2,&n1,&n2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Pf option\n";
			nerror++;
		    }
		    model.addParallel(new
				      FlowLawGBSDry(a1,a2,m1,m2,n1,n2,e1,e2,
						     v1,v2));
		    break;
		}
		case 'g':
                {
                  double a1,a2,m1,m2,n1,n2,r1,r2,e1,e2,v1,v2;
                    if (sscanf(&argv[i][3],
                               "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
                               &a1,&a2,&m1,&m2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 12){
                        cerr << "invalid -Pg option\n";
                        nerror++;
                    }
                    model.addParallel(new
                                      FlowLawGBSWet(a1,a2,m1,m2,n1,n2,r1,r2,e1,e2,
                                                     v1,v2));
                    break;
                }

		default:
		    cerr << "unknown -P option\n";
		    nerror++;
		    break;
		}
		break;

	    case 'S':
	    switch(argv[i][2]){
	    	case 'd':
		{
		  double a1,a2,r1,r2,n1,n2,e1,e2,v1,v2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&r1,&r2,&e1,&e2,&v1,&v2) != 10){
			cerr << "invalid -Sd option\n";
			nerror++;
		    }
		    model.addSequential(new
					FlowLawDisWet(a1,a2,n1,n2,r1,r2,
						      e1,e2,v1,v2));
		    break;
		}
	    	case 'e':
	    	 {
		   double a1,a2,s1,s2,n1,n2,e1,e2;
		    if (sscanf(&argv[i][3],
			       "%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &a1,&a2,&n1,&n2,&s1,&s2,&e1,&e2) != 8){
			cerr << "invalid -Se option\n";
			nerror++;
		    }
		    model.addSequential(new
					FlowLawGen(a1,a2,n1,n2,s1,s2,
						   e1,e2));
		    break;
		}

		default:
		    cerr << "unknown -S option\n";
		    nerror++;
		    break;
		}
		break;

	    case 'F':
	    {
		int i1, i2;
		if (sscanf(&argv[i][2],
			   "%d/%d", &i1, &i2) != 2){
		    cerr << "invalid -F option\n";
		    nerror++;
		}
		ifix.push_back(i1);
		isrc.push_back(i2);
		break;
	    }
	    case 'G':
	    {
	      char* pstr;
	      char s1[MaxStr];
		if (sscanf(&argv[i][2],
			   "%s", s1) != 1){
		    cerr << "invalid -G option\n";
		    nerror++;
		}
		if (!sticky){
		    sgrp.push_back(0);
		    gn++;
		}
		int p;
		pstr = strtok(s1, "/");
		while (pstr!=NULL){
		  gn++;
		  p = atoi(pstr);
		  sgrp.push_back(p);
		  pstr = strtok(NULL, "/");
		}
		sgrp.push_back(0);
		gn++;
		sticky = true;
		break;
	    }
	    case 'B':
	    {
	      new_pairs = true;
	      rfn = &argv[i][2];
	      nNewRuns = model.readNewRunID(rfn);
	      break;
	    }
	    case 'R':
		ran_seed = atoi(&argv[i][2]);
		Parameter::setRanSeed(ran_seed);
		break;
	    case 'C':
	    {
		useCG = true;
		int ntrial = atoi(&argv[i][2]);
		model.useConjugateGradient(ntrial);
		break;
	    }
	    case 'I':
	    {
	      readfile = true;
	      infn = &argv[i][2];
	      FILE *fp;
	      fp = fopen(infn, "r");
	      if (fp==NULL){
		cerr << "Error opening file"; 
		exit (1); 
	      }
	      char output[MaxStr];
	      while(fgets(output, MaxStr, fp)!=NULL){
	      };
	      model.InitModel(readfile,output);
	      sscanf(output, "%6d", &start);
	      break;
	    }
	    case 'i':
	    {
	      // suppliments case 'I' if max_iter is not known, rather no. of iterations that need to be done is given
	      max_iter = max_iter+start;
	    }
	    case 'V':
		verbose = true;
 		break;
	    default:
		cerr << "unknown options\n";
		nerror++;
		break;
	    }
	}else{
	    cerr << "invalid syntax\n";
	    nerror++;
	}
    }

    if (!getD || !getM) nerror++;
    if (model.numParallel()+model.numSequential()==0){

	cerr << "no flow law specified\n";
	nerror++;
    }
    if (nerror){
	cerr << "invalid command option[s] - abort\n";
	exit(1);
    }
    if (verbose){
      cerr << "the number of parallel flow laws: "
	   << model.numParallel() << '\n';
      cerr << "the number of sequential flow laws: "
	   << model.numSequential() << '\n';
    }

    int scalecount = 0;
    for (int i=1; i<=model.numParallel(); i++){
      if (model.parallel(i)->needScaling()){
	scalecount++;
      }
    }
    for(int i=1; i<=model.numSequential(); i++){
      if (model.sequential(i)->needScaling()){
	scalecount++;
      }
    }

    //
    // read experimental data
    //
    model.readData(dfn);
    int ndata = model.numData();
    if (verbose){
      cerr << "the number of data: " << ndata << '\n';
      model.printRunIds(cerr);
    }

    //
    // run MCMC with Gibbs sampling
    //

    long idum1 = long(-(abs(ran_seed)+1)); // for ran2()
    long idum = long(-(abs(ran_seed)+1)); // for ran2()
    model.setUp();
    cerr << "total no. of data-pairs: " << model.rpairs << '\n';
    if (new_pairs){
      cerr << "Number of new data-pairs: " << model.cg_pairs << '\n';
      if (ndata!=nNewRuns)
	error("Size of input new runID and data must be the same\n");
    }

    int np = model.numUnfixedParams();
    Array1d<int> ifixed(np);
    ifixed=0;

    for (int i=1; i<=ifix.size(); i++){
	if (ifix(i)>0 && ifix(i)<=np
	    && isrc(i)>0 && isrc(i)<=np){
	    ifixed(ifix(i)) = isrc(i);
	    if (verbose){
		cerr << model.unfixedParam(ifix(i))->name()
		     << " will be fixed to "
		     << model.unfixedParam(isrc(i))->name()
		     << "\n";
	    }
	}else{
	    error("invalid -F option");
	}
    }

    model.addConstraints(ifixed);
    cerr << "total no. of unfixed parameters (including inter-run biases): ";
    cerr << model.numUnfixedParams() << '\n';

    int nrun = model.runsize();
    int nbias = 0;
    if (nrun>1) nbias = nrun;
    cerr << "no. of inter-run bias: " << nbias << '\n';

    Array1d<double> chi2_vec(max_m*2);
    Array1d<double> calcUnfixedParam; 
    double chi2 = 0;
    double chi2_orig = 0;
    double chi2min = 1e30;
    Array1d<double> prev_minunfixParam; 
    Array1d<int> igrp;
    int gi, pp;

#ifdef _TIMECALC_
    double time_con=0, time_rej=0, time_ran=0, time_CG=0, time_tot=0;
    time_tot = -wtime();
#endif

    for (int ii=start+1; ii<=max_iter; ii++){

	//
	// randomize data
	//

#ifdef _TIMECALC_
      if (ii%dn_ran==0){
	time_ran = 0;
	time_ran = -wtime();
      }
#endif

	int col;
	if (ii%dn_ran==0){
	    col = model.randomizeData();
	    if (outRanData)
	      model.printData(outD);
	}

#ifdef _TIMECALC_
	if (ii%dn_ran==0){
	  time_ran += wtime();
	}
	time_con = 0;
	time_rej = 0;
	time_con = -wtime();
#endif

	//
	// single scan for ordinary model parameters
	//

	int k;
	while (true){
	  k = int(round(ran2(&idum1)*(np-1)+1));
	  if (ifixed(k)>0) k = ifixed(k);
	  if (k<=model.numUnfixedParams()-nbias) break;
	}


	isSticky = false;
	// check if the chosen parameter is a sticky parameter. 
	// If it is, then all sticky parameters grouped with it 
	// will have to be randomized
	if (sticky){
	  int istart, iend;
	  for (int iig=1; iig<=gn; iig++){
	    if (sgrp(iig)==0){
	      if (!isSticky) istart=iig+1;
	      else{
		iend=iig-1;
		break;
	      }
	    }
	    if (sgrp(iig) == k){
	      isSticky = true;
	      //break;
	    }
	  }
	  if (isSticky){
	    gi = iend-istart+1;
	    igrp.resize(gi);
	    for (int iig=istart; iig<=iend; iig++)
	      igrp(iig-istart+1) = sgrp(iig); 
	  }
	}

	// form the prev_min unfixed parameter vector
	if (!isSticky){
	  prev_minunfixParam.resize(3+nbias+scalecount);
	  prev_minunfixParam(1) = model.unfixedParam(k)->value();
	}else{
	  prev_minunfixParam.resize(2+gi+nbias+scalecount);
	  for (int iig=1; iig<=gi; iig++){
	    prev_minunfixParam(iig) = model.unfixedParam(igrp(iig))->value();
	  }
	}

	if (ii-start==1){

	  if (sticky){
	    for (int iig=1; iig<=gi; iig++)
	      model.callBeforeBestFitA(igrp(iig));
	  }
	  chi2 = model.calcChiSq3(k, chi2_orig);
	  chi2min = chi2;
	  
	}

	int p=1;
	if (!isSticky)
	  pp = 1;
	else
	  pp = gi;

	for (p=1; p<=nbias; p++)
	  prev_minunfixParam(pp+p) = model.unfixedParam(model.numUnfixedParams()-nbias+p)->value();
	for (int q=1; q<=model.numParallel(); q++){
	  if (model.parallel(q)->needScaling()){
	    prev_minunfixParam(pp+p) = model.parallel(q)->scaling();
	    p++;
	  }
	}
	for (int q=1; q<=model.numSequential(); q++){
	  if (model.sequential(q)->needScaling()){
	    prev_minunfixParam(pp+p) = model.sequential(q)->scaling();
	    p++;
	  }
        }
	prev_minunfixParam(p+pp) = chi2;
	prev_minunfixParam(p+pp+1) = chi2_orig;

	// If current iteration requires data randomization, randomize data 
	if (ii%dn_ran==0){

	  if (isSticky){
	    for (int iig=1; iig<=gi; iig++){
	      model.callBeforeBestFitA(igrp(iig));
	    }
	  }

	  // recalculate minimum chi2 for randomized data set
	  chi2min = model.calcChiSq3(k, chi2_orig);

	  int p=1;
	  for (p=1; p<=nbias; p++){
	      prev_minunfixParam(pp+p) = model.unfixedParam(model.numUnfixedParams()-nbias+p)->value();
	  }
	  for (int q=1; q<=model.numParallel(); q++){
	    if (model.parallel(q)->needScaling()){
	      prev_minunfixParam(pp+p) = model.parallel(q)->scaling();
	      p++;
	    }
	  }
	  for (int q=1; q<=model.numSequential(); q++){
	    if (model.sequential(q)->needScaling()){
	      prev_minunfixParam(pp+p) = model.sequential(q)->scaling();
	      p++;
	    }
	  }
	  prev_minunfixParam(p+pp) = chi2min;
	  prev_minunfixParam(p+pp+1) = chi2_orig;
	  
	}
	
	// estimate conditional probability distribution
	
	calcUnfixedParam.resize(max_m*(scalecount+nbias+pp));
	
	calcUnfixedParam = 0;
	chi2_vec = 0;
	chi2 = 0;
	chi2_orig = 1e30;
	
	for (int m=1; m<=max_m; m++){

	  if (!isSticky)
	    model.unfixedParam(k)->randomize();
	  else{
	    for (int iig=1; iig<=pp; iig++){
	      model.unfixedParam(igrp(iig))->randomize();
	      model.callBeforeBestFitA(igrp(iig));
	      // Everytime a parameter is randomized, the corresponding 
	      // strain rate prediction (excluding A) has to be recalculated.
	      // When parameter is not sticky, calcChiSq3 function inherently 
	      // performs this function. But, when parameter is sticky, then 
	      // more than 1 parameter is being randomized in the same iteration.
	      // In that case, for each parameter randomized, the corresponding 
	      // flow law willl need to be re-evaluated. callBeforBestFitA() does that.
	    }
	  }
	  	  
	  chi2 = model.calcChiSq3(k,chi2_orig);
	  chi2_vec(2*m-1) = chi2;
	  chi2_vec(2*m) = chi2_orig;
	  
	  if (!isSticky)
	    calcUnfixedParam((nbias+scalecount+pp)*(m-1)+1) = model.unfixedParam(k)->value();
	  else{
	    for (int iig=1; iig<=pp; iig++){
	      calcUnfixedParam((nbias+scalecount+pp)*(m-1)+iig) = model.unfixedParam(igrp(iig))->value();
	    }
	  }
	  
	  int p=1;
	  for (p=1; p<=nbias; p++){
	    calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p) = model.unfixedParam(model.numUnfixedParams()-nbias+p)->value();
	  }
	  for (int q=1; q<=model.numParallel(); q++){
	    if (model.parallel(q)->needScaling()){
	      calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p) = model.parallel(q)->scaling();
	      p++;
	    }
	  }
	  for (int q=1; q<=model.numSequential(); q++){
	    if (model.sequential(q)->needScaling()){
	      calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p) = model.sequential(q)->scaling();
	      p++;
	    }
	  }
	}

#ifdef _TIMECALC_
	time_con += wtime();
	time_rej = -wtime();
	//cerr << "Time taken in Gibb's sampling: " << time_con<<'\n';
#endif 
	
	Array1d<double> sampled_mval_m;
	sampled_mval_m.resize(nbias+scalecount+pp);
	int m_chi2min = 0;
	sampled_mval_m = 0;
	for (int m=1; m<=max_m; m++){ 
	  if (chi2_vec(2*m-1) <= chi2min){
	    chi2min = chi2_vec(2*m-1);
	    m_chi2min = m;
	  }
	}

	if (m_chi2min>0){
	  for (int p=1; p<=(nbias+scalecount+pp); p++)
	    sampled_mval_m(p) = calcUnfixedParam((nbias+scalecount+pp)*(m_chi2min-1)+p); 
	}
	
	// now with the rejection method
	
	int ir=0;
	while (true){
	  
	  int m = int(round(ran2(&idum)*(max_m-1)+1));
	  ir++;
	  
	  double prob = exp(-0.5*(chi2_vec(2*m-1)-chi2min));
	  if (ran2(&idum)<prob){
	    chi2 = chi2_vec(2*m-1);
	    chi2_orig = chi2_vec(2*m);
	    if (!isSticky)
	      model.unfixedParam(k)->setValue(calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp));
	    else{
	      for (int iig=1; iig<=pp; iig++)
		model.unfixedParam(igrp(iig))->setValue(calcUnfixedParam((nbias+scalecount+pp)*(m-1)+iig));
	    }
	    int p=1;
	    for (p=1; p<=nbias; p++){
	      model.unfixedParam(model.numUnfixedParams()-nbias+p)->setValue(calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p));
	    }
	    for (int q=1; q<=model.numParallel(); q++){
	      if (model.parallel(q)->needScaling()){
		model.parallel(q)->setScaling(calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p));
		p++;
	      }
	    }
	    for (int q=1; q<=model.numSequential(); q++){
	      if (model.sequential(q)->needScaling()){
		model.sequential(q)->setScaling(calcUnfixedParam((nbias+scalecount+pp)*(m-1)+pp+p));
		p++;
	      }
	    }	     
	    break;
	  }
	  if (ir>max_m){ // this extra if-clause is for efficiency
	    
	    if (m_chi2min>0){
	      
	      //cerr << "acceptance cond. not met, but min. model exists\n";
	      
	      chi2 = chi2_vec(2*m_chi2min-1);
	      chi2_orig = chi2_vec(2*m_chi2min);
	      if (!isSticky)
		model.unfixedParam(k)->setValue(sampled_mval_m(1));
	      else{
		for (int iig=1; iig<=pp; iig++)
		  model.unfixedParam(igrp(iig))->setValue(sampled_mval_m(iig));
	      }
	      
	      int p=1;
	      for (p=1; p<=nbias; p++)
		model.unfixedParam(model.numUnfixedParams()-nbias+p)->setValue(sampled_mval_m(pp+p));
	      for(int q=1; q<=model.numParallel(); q++){
		if (model.parallel(q)->needScaling()){
		  model.parallel(q)->setScaling(sampled_mval_m(pp+p));
		  p++;
		}
	      }
	      for (int q=1; q<=model.numSequential(); q++){
		if (model.sequential(q)->needScaling()){
		  model.sequential(q)->setScaling(sampled_mval_m(pp+p));
		  p++;
		}
	      }
	      break;
	      
	    }else{
	      //cerr << "Acceptance condition not met and no min. model.";
	      //cerr << " chi2 of prev model:" << prev_minunfixParam(3) << "\n";
	      if (!isSticky)
		model.unfixedParam(k)->setValue(prev_minunfixParam(1));
	      else{
		for (int iig=1; iig<=pp; iig++)
		  model.unfixedParam(igrp(iig))->setValue(prev_minunfixParam(iig));
	      }
	      
	      int p=1;
	      for (p=1; p<=nbias; p++){
		model.unfixedParam(model.numUnfixedParams()-nbias+p)->setValue(prev_minunfixParam(pp+p));
	      }
	      for(int q=1; q<=model.numParallel(); q++){
		if (model.parallel(q)->needScaling()){
		  model.parallel(q)->setScaling(prev_minunfixParam(pp+p));
		  p++;
		}
	      }
	      for (int q=1; q<=model.numSequential(); q++){
		if (model.sequential(q)->needScaling()){
		  model.sequential(q)->setScaling(prev_minunfixParam(pp+p));
		  p++;
		}
	      }
	      chi2 = prev_minunfixParam(p+pp);
	      chi2_orig = prev_minunfixParam(p+pp+1);
	      break;
	    }
	  }
	}
	
#ifdef _TIMECALC_
	time_rej += wtime();
	//cerr << "Time taken in rejection method: " << time_rej << " ";
#endif	
	
	Array1d<double> setParams(model.numUnfixedParams()+scalecount);
	for (p=1; p<=model.numUnfixedParams(); p++){
	  setParams(p) = model.unfixedParam(p)->value();
	}
	for(int q=1; q<=model.numParallel(); q++){
	  if (model.parallel(q)->needScaling()){
	    setParams(p) = model.parallel(q)->scaling();
	    p++;
	  }
	}
	for (int q=1; q<=model.numSequential(); q++){
	  if (model.sequential(q)->needScaling()){
	    setParams(p) = model.sequential(q)->scaling();
	    p++;
	  }
	}
	for (p=1; p<=model.numUnfixedParams(); p++)
	  model.unfixedParam(p)->setValue(setParams(p));
	
	for(int q=1; q<=model.numParallel(); q++){
	  if (model.parallel(q)->needScaling()){
	    model.parallel(q)->setScaling(setParams(p));
	    p++;
	  }
	}
	for (int q=1; q<=model.numSequential(); q++){
	  if (model.sequential(q)->needScaling()){
	    model.sequential(q)->setScaling(setParams(p));
	    p++;
	  }
	}
	model.set_fmat();
	
	// output the current solution
	if (ii%dn_out == 0){
	  char line[MaxStr];
	  sprintf(line, "%6d %6.5e %6.5e ", ii, chi2/(model.rpairs), chi2_orig/(model.rpairs));
	  cout << line;

	  for (int i=1; i<=model.numUnfixedParams(); i++){
	    sprintf(line, "%6.5e ", model.unfixedParam(i)->value());
	    cout << line;
	  }
	  
	  for (int i=1; i<=model.numParallel(); i++){
	    if (model.parallel(i)->needScaling()){
	      sprintf(line, "%6.5e ", model.parallel(i)->scaling());
	      cout << line;
	    }
	  }
	  
	  for (int i=1; i<=model.numSequential(); i++){
	    if (model.sequential(i)->needScaling()){
	      sprintf(line, "%6.5e ", model.sequential(i)->scaling());
	      cout << line;
	    }
	  }
	  cout << '\n';
	  cout.flush();
	}	
    }
#ifdef _TIMECALC_
    time_tot += wtime();
    cerr << "total time: " << time_tot << "s = " << time_tot/3600 << "hrs.\n";
#endif
}
