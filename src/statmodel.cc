/*
 * statmodel_stickyPair2.cc
 *
 * Jun Korenaga
 * Summer 2008
 *
 * CJ: 
 * cost function edited to remove log - 23 July 2015
 * September 2015 : fit model to log(e2) - log(e1) so that there is no need for * inter-run bias
 * December 2017 : Testing reduced number of pairs for CG search.
 *                 Also, incorporating a solution for sticky parameters. 
 * February 2018 : Corrected code to incorporate flow laws with fixed A.
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <ctime>
#include "labmc_nr.h"
#include "statmodel.h"
#include "util.h"

//#define _CHECK_

const double StatModel::max_chi2 = 1e30;

StatModel::StatModel()
{
    do_CG = false;
    is_constrained = false;
    infile = false;
    less_pairs = false;
}

StatModel::~StatModel()
{
    for (int i=1; i<=para.size(); i++){
        delete para(i);
	para(i) = nullptr;
    }
    for (int i=1; i<=seq.size(); i++){
        delete seq(i);
	seq(i) = nullptr;
    }
}

void StatModel::readData(char *fn)
{
    int ndata = countLines(fn);
    orig_data.resize(ndata);
    data.resize(ndata);
    runID.resize(ndata);

    istream* pin = new ifstream(fn);
    for (int i=1; i<=ndata; i++){
        double T, dT, p, dp, e, de, sig, dsig;
        double id;
        *pin >> T >> dT >> p >> dp >> e >> de >> sig >> dsig;
        orig_data(i).setTemperature(T,dT);
        orig_data(i).setPressure(p*1e9,dp*1e9);
        orig_data(i).setStrainRate(e,de);
        orig_data(i).setStress(sig,dsig);

        for (int j=1; j<=more_state.size(); j++){
            double val,dval;
            *pin >> val >> dval;
            setMoreState(orig_data(i),more_state(j),val,dval);
        }

        *pin >> id;
        runID(i) = int(id);
    }

    data = orig_data;
}

int StatModel::randomizeData()
{
    for (int i=1; i<=orig_data.size(); i++){
        double T = orig_data(i).temperature();
        double dT = orig_data(i).temperatureError();
        double p = orig_data(i).pressure();
        double dp = orig_data(i).pressureError();
        double e = orig_data(i).strainRate();
        double de = orig_data(i).strainRateError();
        double sig = orig_data(i).stress();
        double dsig = orig_data(i).stressError();

        double newT = randomizeValue(T,dT);
        double newp = randomizeValue(p,dp);
        double newe = randomizeValue(e,de);
        double newsig = randomizeValue(sig,dsig);
        
        data(i).setTemperature(newT,dT);
        data(i).setPressure(newp,dp);
        data(i).setStrainRate(newe,de);
        data(i).setStress(newsig,dsig);

        double val, dval, newval;
        for (int j=1; j<=more_state.size(); j++){
            switch(more_state(j)){
            case State::grain_size:
            val = orig_data(i).grainSize();
            dval = orig_data(i).grainSizeError();
            newval = randomizeValue(val,dval);
            data(i).setGrainSize(newval,dval);
            break;
            case State::water_content:
            val = orig_data(i).waterContent();
            dval = orig_data(i).waterContentError();
            newval = randomizeValue(val,dval);
            data(i).setWaterContent(newval,dval);
            break;
            case State::oxygen_fugacity:
            val = orig_data(i).oxygenFugacity();
            dval = orig_data(i).oxygenFugacityError();
            newval = randomizeValue(val,dval);
            data(i).setOxygenFugacity(newval,dval);
            break;
            default:
            error("StatModel::randomizeData - invalid istate");
            break;
            }
        }
    }

    beforeBestFitA(-1); // recalculate with newly randomized data
    int col = 8+2*more_state.size();
    return col;
}

double StatModel::randomizeValue(double v, double dv)
{
//    double newv = v + dv*gasdev(&data_idum);
    double newv = v + dv*2.0*(ran2(&data_idum)-0.5);
    //    cerr << " r_params: " << v << " " << dv << " " << newv << '\n';
    return newv;
}

Array1d<double> StatModel::outData()
{
  Array1d<double>ranData(data.size()*(8+2*more_state.size()));
  int r=data.size();
  
  for (int i=1; i<=data.size(); i++){
    ranData(i) = data(i).temperature();
    ranData(r+i) = data(i).temperatureError();
    ranData(2*r+i) = data(i).pressure();
    ranData(3*r+i) = data(i).pressureError();
    ranData(4*r+i) = data(i).strainRate();
    ranData(5*r+i) = data(i).strainRateError();
    ranData(6*r+i) = data(i).stress();
    ranData(7*r+i) = data(i).stressError();
    
    for (int j=1; j<=more_state.size(); j++){
    switch(more_state(j)){
    case State::grain_size:
      ranData((6+2*j)*r+i) = data(i).grainSize();
      ranData((7+2*j)*r+i) = data(i).grainSizeError();
      break;
    case State::water_content:
      ranData((6+2*j)*r+i) = data(i).waterContent();
      ranData((7+2*j)*r+i) = data(i).waterContentError();
      break;
    case State::oxygen_fugacity:
      ranData((6+2*j)*r+i) = data(i).oxygenFugacity();
      ranData((7+2*j)*r+i) = data(i).oxygenFugacityError();
      break;
    default:
      error("StatModel::outData - invalid istate");
      break;
    }
    }
  }
  return ranData;
}

void StatModel::copyData(Array1d<double>& ranData)
{
  int r = data.size();
    for (int i=1; i<=data.size(); i++){
      double T = ranData(i);
      double dT = ranData(r+i);
      double p = ranData(2*r+i);
      double dp = ranData(3*r+i);
      double e = ranData(4*r+i);
      double de = ranData(5*r+i);
      double sig = ranData(6*r+i);
      double dsig = ranData(7*r+i);

      data(i).setTemperature(T,dT);
      data(i).setPressure(p,dp);
      data(i).setStrainRate(e,de);
      data(i).setStress(sig,dsig);

      double val, dval;
      for (int j=1; j<=more_state.size(); j++){
	switch(more_state(j)){
	case State::grain_size:
	  val = ranData((6+2*j)*r+i);
	  dval = ranData((7+2*j)*r+i);
	  data(i).setGrainSize(val,dval);
	  break;
	case State::water_content:
	  val = ranData((6+2*j)*r+i);
	  dval = ranData((7+2*j)*r+i);
	  data(i).setWaterContent(val,dval);
	  break;
	case State::oxygen_fugacity:
	  val = ranData((6+2*j)*r+i);
	  dval = ranData((7+2*j)*r+i);
	  data(i).setOxygenFugacity(val,dval);
	  break;
	default:
	  error("StatModel::copyData - invalid istate");
	  break;
            }
      }
    }
    beforeBestFitA(-1); // recalculate with newly randomized data
}

void StatModel::printData(char *outD)
{
  //  odf << "data:\n";
  std::ofstream odf (outD, std::ofstream::app);

  for (int i=1; i<=data.size(); i++){ //i+=10){
    odf << i << " " << data(i).temperature() << " " << data(i).temperatureError() << " ";
    odf << data(i).pressure()*1e-9 << " " << data(i).pressureError()*1e-9 << " ";
    odf << data(i).strainRate() << " " << data(i).strainRateError() << " ";
    odf << data(i).stress() << " " << data(i).stressError();
        for (int j=1; j<=more_state.size(); j++){
	    odf << " ";
	    switch(more_state(j)){
            case State::grain_size:
	      odf << data(i).grainSize() << " " << data(i).grainSizeError() << " ";
            break;
            case State::water_content:
	      odf << data(i).waterContent() << " " << data(i).waterContentError() << " ";
            break;
            case State::oxygen_fugacity:
	      odf << data(i).oxygenFugacity() << " " << data(i).oxygenFugacityError() << " ";
            break;
            default:
            error("StatModel::randomizeData - invalid istate");
            break;
            }
        }
	odf << runID(i) << '\n';
  }
  //  os << '\n';
}

void StatModel::printRunIds(ostream& os)
{
    os << "runids: ";
    for (int i=1; i<=runID.size(); i++){
	os << runID(i) << " ";
    }
    os << '\n';
}

int StatModel::readNewRunID(char *rfn)
{
    less_pairs = true;
    int ndata = countLines(rfn);
    newID.resize(ndata);

    istream* pin = new ifstream(rfn);
    for (int i=1; i<=ndata; i++){
         double id;
        *pin >> id;
        newID(i) = int(id);
    }
    return newID.size();
}

void StatModel::stateToRead(int i)
{
    more_state.push_back(i);
}

void StatModel::setMoreState(State& d, int istate, 
			     double val, double dval)
{
    switch(istate){
    case State::grain_size:
	d.setGrainSize(val,dval);
	break;
    case State::water_content:
	d.setWaterContent(val,dval);
	break;
    case State::oxygen_fugacity:
	d.setOxygenFugacity(val,dval);
	break;
    default:
	error("StatModel::setMoreState - invalid istate");
	break;
    }
}

/*void StatModel::setupBiasCorrection(bool set_bias, double bmin, double bmax)
//void StatModel::setupBiasCorrection()
{
  //    do_bias = true;
  do_bias = set_bias;
      min_bias = bmin;
      max_bias = bmax;
}
*/

void StatModel::setUp()
{
    //
    // check if all of required states are given
    //
    for (int i=1; i<=para.size(); i++){
	if (para(i)->checkStates(more_state)==false){
	    error("StatModel::setUp - missing state(s)");
	}
    }
    for (int i=1; i<=seq.size(); i++){
	if (seq(i)->checkStates(more_state)==false){
	    error("StatModel::setUp - missing state(s)");
	}
    }

    //
    // classify flow laws
    //
    i_para_scaling.resize(0);
    i_para_noscaling.resize(0);
    i_seq_scaling.resize(0);
    i_seq_noscaling.resize(0);

    for (int i=1; i<=para.size(); i++){
	if (para(i)->needScaling()){
	    i_para_scaling.push_back(i);
	}else{
	    i_para_noscaling.push_back(i);
	}
    }
    for (int i=1; i<=seq.size(); i++){
	if (seq(i)->needScaling()){
	    i_seq_scaling.push_back(i);
	}else{
	    i_seq_noscaling.push_back(i);
	}
    }

    // allocate relevant arrays
    fmat_para.resize(orig_data.size(),para.size());
    fmat_seq.resize(orig_data.size(),seq.size());
    edot_para.resize(orig_data.size(),para.size());
    edot_seq.resize(orig_data.size(),seq.size());
    edot_para_tmp.resize(orig_data.size(),para.size());
    edot_seq_tmp.resize(orig_data.size(),seq.size());
    AtA1.resize(i_para_scaling.size(),i_para_scaling.size());
    AtA2.resize(i_seq_scaling.size(),i_seq_scaling.size());
    Atb1.resize(i_para_scaling.size(),1);
    Atb2.resize(i_seq_scaling.size(),1);
    chivec.resize(orig_data.size());
    int n_scaling=i_para_scaling.size()+i_seq_scaling.size();
    logA.resize(n_scaling);
    maxlogA.resize(n_scaling);
    minlogA.resize(n_scaling);
    bestlogA.resize(n_scaling);
    logBbyA.resize(n_scaling);
    maxlogBbyA.resize(n_scaling);
    minlogBbyA.resize(n_scaling);
    bestlogBbyA.resize(n_scaling);
    BbyA.resize(n_scaling);
    grad.resize(n_scaling);
    new_grad.resize(n_scaling);
    direc.resize(n_scaling);
    tmp_logBbyA.resize(n_scaling);
    tmp_logA.resize(n_scaling);
    //
    // count unfixed parameters
    //
    unfixed_params.resize(0);
    param_type.resize(0);
    param_id.resize(0);

    for (int i=1; i<=para.size(); i++){
	for (int j=1; j<=para(i)->numUnfixedParams(); j++){
	  unfixed_params.push_back(para(i)->unfixedParam(j));
	  param_type.push_back(type_para);
	  param_id.push_back(i);
	}
    }
    for (int i=1; i<=seq.size(); i++){
        for (int j=1; j<=seq(i)->numUnfixedParams(); j++){
	  unfixed_params.push_back(seq(i)->unfixedParam(j));
	  param_type.push_back(type_seq);
	  param_id.push_back(i);
	}
    }

    Array1d<int> tmpid(runID.size());
    tmpid = runID;
    //tmpid.sort_unique(); 
    // sort_unique picks the unique run IDs & sorts them in increasing order
    tmpid.f_unique(); // this just picks the unique run IDs; doesn't sort
    for (int i=1; i<=tmpid.size(); i++){
      uid[tmpid(i)] = i;
    }
    nbias = tmpid.size();
    if (nbias == 1) bias.resize(nbias-1);
    else bias.resize(nbias);     

    if (nbias>1){
      for (int i=1; i<=nbias; i++){
	char str[MaxStr];
	sprintf(str,"Bias[%d]",i);
	bias(i) = new Parameter(0.0,str);
	unfixed_params.push_back(bias(i));
	param_type.push_back(type_bias);
	param_id.push_back(i);
      }
    }

    int n;
    rpairs = 0;
    for (int i=1; i<=nbias; i++){
	  n = 0;
	  for (int j=1; j<=runID.size(); j++){
	    if (runID(j)==tmpid(i)){
	      n++;
	    }
	  }
	  rpairs += n*(n-1)/2;
    }

    // If fewer pairs have to be made in CG search, then count those.
   Array1d<int> ntmpid(newID.size());
    ntmpid = newID;
    //    ntmpid.sort_unique();
    ntmpid.f_unique();
    for (int i=1; i<=ntmpid.size(); i++){
      uid[ntmpid(i)] =  i;
    }

    if (do_CG && less_pairs){
      int j, m, pri, ii=1, jj=1;
      cg_pairs = 0;
      for (int i=1; i<=nbias; i++){
	m = 0;
	pri = ii;
	for (j=jj; j<=data.size(); j++){
	  if (runID(j)==tmpid(i)){
	    if (newID(j)==ntmpid(ii))
	      m++;
	    else{
	      cg_pairs += m*(m-1)/2;
	      ii++;
	      m=1;
	    }
	  }else{
	    cg_pairs += (m*(m-1)/2)+(2*(ii-pri+1)*(ii-pri));
	    ii++;
	    break;
	  }
	}
	jj = j;
      }
      cg_pairs += m*(m-1)/2 + 2*(ii-pri+1)*(ii-pri);
    }

    // If this inversion is a continuation of a previous output, 
    // read the last model from the old output file.
    if (infile){
      char * instr;
      int nparams;
      if (nbias>1){
        nparams = numUnfixedParams() - nbias;
      }else{
	nparams = numUnfixedParams();
      }

      int tok=1;
      instr = strtok(inLine, " ");
      while (instr!=NULL){
	if (tok==3+numUnfixedParams()) break;
	tok = tok+1;
	instr  = strtok(NULL, " ");
	if (tok>3){
	  double inval = atof(instr);
	  if (tok<=3+nparams){
	    unfixed_params(tok-3)->setValue(inval);
	  }else{
	    bias(tok-3-nparams)->setValue(inval);
	  }
	}
      }
    }

    // calculate size of Fmat vectors according to the number
    // the number of data pairs to be evaluated within a run
    int indx = 1;
    int rundatamax = 0; 
    for (int i=1; i<=nbias; i++){
      int k = 0;
      for (int j=indx; j<=orig_data.size(); j++){
	if (runID(j)==runID(indx)) k++;
	else{
	  if (rundatamax < k) rundatamax = k;
	  indx = j;
	  break;
	}
      }
    }
    Fmat1.resize(rundatamax,i_para_scaling.size());
    Fmat2.resize(rundatamax,i_seq_scaling.size());
    tmp_logAX.resize(nbias);
    //
    // data-related pre-calculation
    //
    inv_dedot.resize(orig_data.size());
    bvec.resize(orig_data.size());
    tmp_bvec.resize(orig_data.size());
    data_lognorm2=0.0;
    
    //cerr << "\nbvec initialization:\n";

    for (int i=1; i<=orig_data.size(); i++){
        double val1 = data(i).strainRateError();
        inv_dedot(i) = 1.0/val1;
        bvec(i) = data(i).strainRate();
	//cerr << bvec(i) << '\n';
	double val2 = log(data(i).strainRate());  
	data_lognorm2 += val2*val2; // [log(edot)]^2
    }

    isrc.resize(unfixed_params.size());
    for (int i=1; i<=isrc.size(); i++){
	isrc(i).resize(0);
	isrc(i).push_back(i);
    }
    //
    // initialize chi2 calculation
    //
    beforeBestFitA(-1); // calculates  edot_predicted/pre-exponential_factor 

    cg_idum = long(-(abs(data_lognorm2*1000)+1)); // for ran2();
    data_idum = long(-(abs(data_lognorm2*1000)+2)); // for ran2();
}

void StatModel::addConstraints(const Array1d<int>& ifixed)
{
    for (int i=1; i<=ifixed.size(); i++){
	if (ifixed(i)>0) isrc(ifixed(i)).push_back(i);
    }
    is_constrained=true;
    fixParams();

}

void StatModel::fixParams()
{
  if (is_constrained){
    for (int i=1; i<=isrc.size(); i++){
       if (isrc(i).size()>1){
	 for (int j=2; j<=isrc(i).size(); j++){
	    unfixed_params(isrc(i)(j))->setValue(unfixed_params(isrc(i)(1))->value());
	  }
       }
    }
  }
}

void StatModel::beforeBestFitA(int k)
{
    if (k<0){
	for (int i=1; i<=para.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		fmat_para(j,i) = para(i)->prepPredict(data(j));
	    }      
	}

	for (int i=1; i<=seq.size(); i++){
	    for (int j=1; j<=data.size(); j++){
		fmat_seq(j,i) = seq(i)->prepPredict(data(j));
	    }      
	}
    }else{
	for (int i=1; i<=isrc(k).size(); i++){
	    int kk = isrc(k)(i);

	    switch(param_type(kk)){
	    case type_para:
	    {
		int i = param_id(kk);
		for (int j=1; j<=data.size(); j++){
		    fmat_para(j,i) = para(i)->prepPredict(data(j));
		}      
		break;
	    }
	    case type_seq:
	    {
		int i = param_id(kk);
		for (int j=1; j<=data.size(); j++){
		    fmat_seq(j,i) = seq(i)->prepPredict(data(j));
		}      
		break;
	    }
	    case type_bias:
	      //*// do nothing
		break;
	    }
	}
    }
}

bool StatModel::calcBestFitA3()
{
  //*// 1. para_scaling==0 && seq_scaling==0 
    //    (this is not an impossible case. other parameters can vary)
    // 2. para_scaling>0 && seq_scaling>0 
    //    ---> do_CG
    // 3. para_scaling>0 && seq_scaling==0
    //    ---- but seq_nosclaing can still be non-zero
    // 4. para_scaling==0 && seq_scaling>0 
    //    ---- but para_noscaling can still be non-zero
  //*

    if (i_para_scaling.size()==0 && i_seq_scaling.size()==0){
      //*// need to do nothing
        return true;
    }
    
     if (do_CG) {
        if (i_para_scaling.size()>0 && i_seq_scaling.size()>0){
	  if (i_para_noscaling.size()>0 || i_seq_noscaling.size()>0) 
	    return calcBestFitA_CG2();
	  else return calcBestFitA_CG3();
        }else{
	  if (i_para_scaling.size()>1 || i_seq_scaling.size()>1){
            if (i_para_noscaling.size()>0 || i_seq_noscaling.size()>0) 
	      return calcBestFitA_CG2();
	    else return calcBestFitA_CG3();
	  }else{
            cerr << " CG cannot invert for only 1 flow law\n";
            exit(1);
	  }
	}
    }else{
       if (i_para_scaling.size()>1 || i_seq_scaling.size()>1){
	 cerr << "Enable do_CG is recommended\n";
	 exit(1);
       }
     }
    
     if (i_para_scaling.size()>0 && i_seq_scaling.size()>0){
       if (!do_CG){
	 cerr << "do_CG needs to be set\n";
	 exit(1);
       }
     }

    if (i_para_scaling.size()>0 && i_seq_scaling.size()==0){
       // i_seq_scaling.size() == 0 but
       // i_seq_noscaling.size() may not be zero.
       //
       // edot = x1*f1+x2*f2+...+e1+e2+... + 1/(1/g1+1/g2+...)
       // --> edot-(e1+e2+...)-(1/(1/g1+..)) = x1*f1+x2*f2+...
       // --> tmp_bvec = Fmat*xvec
       //

       // prepare tmp_bvec
       tmp_bvec = bvec;
       if (i_para_noscaling.size()>0){
	  for (int i=1; i<=i_para_noscaling.size(); i++){
	    double Afix = para(i_para_noscaling(i))->scaling();
             for (int j=1; j<=data.size(); j++){
	       tmp_bvec(j) -= fmat_para(j,i_para_noscaling(i))*Afix;
	     }
	  }
       }

       if (i_seq_noscaling.size()>0){
	  double val, val1;
	  for (int j=1; j<=data.size(); j++){
             val = 0.0;
             for (int i=1; i<=i_seq_noscaling.size(); i++){
	       double Afix = seq(i_seq_noscaling(i))->scaling();
	       val += 1./(fmat_seq(j,i_seq_noscaling(i))*Afix);
             }
	     tmp_bvec(j) -= 1/val;
          }
       }


       if (removeNonPositive(tmp_bvec)==false) return false;
       Array1d<double> tmp_bvec2 = tmp_bvec;
 
       if (i_para_scaling.size()>1){
	 // divide tmp_bvec with error in strain rate
	 for (int j=1; j<=data.size(); j++){
	   tmp_bvec(j) = tmp_bvec(j)*inv_dedot(j);
	 } 

          // preprare Fmat
	 int indx = 1, i, j, k, l;
	  Array1d<double> logAtmp(i_para_scaling.size());

	  for (i=1; i<=i_para_scaling.size(); i++){
	    double diff = 0.0;
	    for (j=1; j<=data.size(); j++){
	      diff = log(bvec(j)) - log(fmat_para(j,i_para_scaling(i)));
	      if (j==1) logAtmp(i) = diff;
	      else{
		if (logAtmp(i) > diff) logAtmp(i) = diff;
	      }
	    }
	  }

	  int FLAG = 1;
	  while (FLAG){
	    k = 0;
      	    for (j=indx; j<=data.size(); j++){
	       if (runID(j)==runID(indx)){
		  k++;
                  for (i=1; i<=i_para_scaling.size(); i++){
                     Fmat1(k,i) = fmat_para(j,i_para_scaling(i))*inv_dedot(j)*exp(logAtmp(i));
                  }
               }else break;
            }
            int flen = k;

            // solve linear system
	    for (k=1; k<=i_para_scaling.size(); k++){
              for (l=1; l<=k; l++){
                 AtA1(k,l) = 0.0;
                 for (i=1; i<=flen; i++){
                   AtA1(k,l) += Fmat1(i,k)*Fmat1(i,l);
                 }
                 AtA1(l,k) = AtA1(k,l);
              }
              Atb1(k,1) = 0.0;
              for (i=1; i<=flen; i++){
                Atb1(k,1) += Fmat1(i,k)*tmp_bvec(j-1-flen+i);
              }
	    }
        
            gaussj(AtA1.toRecipe(), i_para_scaling.size(), Atb1.toRecipe(), 1);

	    // check if there's a negative scaling constant...
	    for (i=1; i<=i_para_scaling.size(); i++){
	      if (Atb1(i,1)<0) {
		return false;
	      }else {
		double sumBbyA = exp(logAtmp(i)-logAtmp(1))*Atb1(i,1)/Atb1(1,1); 
		if (indx==1) BbyA(i) = sumBbyA/nbias;
		else BbyA(i) += sumBbyA/nbias;
	      }
	    }

	    indx = j;
	    if (j>=data.size()) FLAG = 0;
	  }
	  /* Atb1 contains the scaling constants A & B... obtained 
	    from LSQR fitting for given run
            We calculate B/A for each run. logBbyA contains the 
	    average log(BbyA) */
	  AXfromBbyA(tmp_bvec2);
       }else{ 
	// i.e. case where i_para_scaling.size()==1 && i_seq_scaling.size()==0
	 BbyA(1) = 1;
	 AXfromBbyA(tmp_bvec2);
       }

       for (int i=1; i<=i_para_scaling.size(); i++){
	 para(i_para_scaling(i))->setScaling(exp(logA(i)));
	 //cerr << "logA:" << logA(1) << " :" << para(i_para_scaling(i))->scaling() << "\n";
       }
    }

    if (i_seq_scaling.size()>0 && i_para_scaling.size()==0){
	// but i_para_noscaling.size() may not be zero.
	//
	// edot = e1+e2+.. 1/(1/g1+1/g2+... +1/x1*h1+1/x2*h2+...)
	// 1/(edot-(e1+e2+...)) - (1/g1+1/g2+...) = 1/x1*h1+1/x2*h2+...
	//                                        = y1/h1+y2/h2+...
	// ---> tmp_bvec = Fmat*yvec
	//

       // prepare tmp_bvec
       tmp_bvec = bvec;
       if (i_para_noscaling.size()>0){
        for (int i=1; i<=i_para_noscaling.size(); i++){
	  double Afix = para(i_para_noscaling(i))->scaling();
          for (int j=1; j<=data.size(); j++){
	    tmp_bvec(j) -= fmat_para(j,i_para_noscaling(i))*Afix;
          }
        }
       }
       if (removeNonPositive(tmp_bvec)==false) return false;

       for (int j=1; j<=data.size(); j++){
          tmp_bvec(j) = 1.0/tmp_bvec(j);
       }

       if (i_seq_noscaling.size()>0){
	 for (int j=1; j<=data.size(); j++){
          for (int i=1; j<=i_seq_noscaling.size(); i++){
	    double Afix = seq(i_seq_noscaling(i))->scaling();
	    tmp_bvec(j) -= 1.0/(fmat_seq(j,i_seq_noscaling(i))*Afix);
          }
        }
       }
       if (removeNonPositive(tmp_bvec)==false) return false;

       Array1d<double> tmp_bvec2 = tmp_bvec;

       // define tmp_bvec
       for (int j=1; j<=data.size(); j++){
	 tmp_bvec(j) = tmp_bvec(j)*data(j).strainRateError();
       }

	// prepare Fmat
       int indx = 1;
       int i, j, k, l;
       int FLAG = 1;

       while (FLAG){
	 k=0;
	 for (j=indx; j<=data.size(); j++){
	     if (runID(j)==runID(indx)){
	       k++;
	       for (i=1; i<=i_seq_scaling.size(); i++){
	 	  Fmat2(k,i) =
		     data(j).strainRateError()/fmat_seq(j,i_seq_scaling(i));
	       }
	     }else break;
	 }

	 int flen = k;

	 // solve linear system
	 for (k=1; k<=i_seq_scaling.size(); k++){
	    for (l=1; l<=k; l++){
		AtA2(k,l) = 0.0;
		for (i=1; i<=flen; i++){
		    AtA2(k,l) += Fmat2(i,k)*Fmat2(i,l);
		}
		AtA2(l,k) = AtA2(k,l);
	    }
	    Atb2(k,1) = 0.0;
	    for (i=1; i<=flen; i++){
		Atb2(k,1) += Fmat2(i,k)*tmp_bvec(j-1-flen+i);
	    }
	 }

	 gaussj(AtA2.toRecipe(), i_seq_scaling.size(), Atb2.toRecipe(), 1);

	 // check if there's a negative scaling constant...
	 for (i=1; i<=i_seq_scaling.size(); i++){
	   if (Atb2(i,1)<0){
	     return false;
	   }else{
	     /* double sumBbyA = Atb2(i,1)/Atb2(1,1); // I made this change in Feb, 2018
	      if (indx==1) BbyA(i) = sumBbyA/nbias; // Not sure if it is correct though.
	      else BbyA(i) += sumBbyA/nbias;*/
	     double sumBbyA = Atb2(1,1)/Atb2(i,1);
	     if (indx==1) BbyA(i) = sumBbyA/nbias;
	     else BbyA(i) += sumBbyA/nbias;
	    }
	 }
	 indx = j;
	 if (j>=data.size()) FLAG=0;
       }

       AXfromBbyA(tmp_bvec2);

       for (int i=1; i<=i_seq_scaling.size(); i++){
	 seq(i_seq_scaling(i))->setScaling(exp(logA(i)));
       }
    }
    
    //    cerr << "scaling coeff: " << para(i_para_scaling(1))->scaling();
    return true;
}

void StatModel::AXfromBbyA(Array1d<double>& tmp_bvec_)
{
  //  afterBestBbyA();

  //cerr << "edot  ediff  bvec  egbs/A   log(bvec)-log(egbs/A)\n";

  int indx = 1;
  int j, u=0;
  double diff1, diff2=0.0;
  int FLAG = 1;
  while (FLAG){
    diff1 = 0.0;
    int k=0;
    for (j=indx; j<=data.size(); j++){
      if (runID(j)==runID(indx)){
	 k++;
	 double edot_total_tmp = 0.0;

	 if (i_para_scaling.size()>0 && i_seq_scaling.size()==0){
	   for (int i=1; i<=i_para_scaling.size(); i++){
	     edot_total_tmp += fmat_para(j,i_para_scaling(i))*BbyA(i);
	   }
	 }
	 
	 if (i_para_scaling.size()==0 && i_seq_scaling.size()>0){
	   // this is in general incorrect
	   double val = 0;
	   for (int i=1; i<=i_seq_scaling.size(); i++){
	     edot_total_tmp += 1.0/(fmat_seq(j,i_seq_scaling(i))*BbyA(i));
	     // this is actually inverse of edot, but tmp_bvec_ in this case is 
	     // defined as the inverse
	   }
	 }

	 diff1 += log(tmp_bvec_(j)) - log(edot_total_tmp);

	 /*cerr << data(j).strainRate() << " ";
	 cerr << fmat_para(j,i_para_noscaling(1))*para(i_para_noscaling(1))->scaling() << " ";
	 cerr << tmp_bvec_(j) << " " << edot_total_tmp << " ";
	 cerr << log(tmp_bvec_(j)) - log(edot_total_tmp) << '\n';*/
      
      }else break; 
    }
    
    if (i_para_scaling.size()==0 && i_seq_scaling.size()>0)
      diff1 *= -1;
    // log(1/edot) - log(1/fmat_seq) = log(1/Aexp(X)) = -X-log(A)
    
    if (nbias>1){
      u++;
      tmp_logAX(u) =  diff1/k;
    }
    
    diff2 += diff1;
    indx = j;
    if (j>=data.size()) FLAG = 0;
  }

  //  cerr << "diff2:" << diff2 << "diff2/N=" << diff2/data.size() << '\n';

  for (int ii=1; ii<=BbyA.size(); ii++){
    logA(ii) = log(BbyA(ii)) + (diff2/data.size());
    //cerr << "logA:" << logA(ii) << " ";
  }
  
  if (nbias>1){
    for (int u=1; u<=nbias; u++){
      bias(u) -> setValue(tmp_logAX(u)-(diff2/data.size()));
    }
  }

}

bool StatModel::removeNonPositive(Array1d<double>& v)
{
    double vmin=1e30;
    for (int j=1; j<=v.size(); j++){
      if (v(j)>0 && v(j)<vmin){
        vmin = v(j);
      }
    }

    if (vmin==1e30){ // i.e., all tmp_bvec is negative
      return false;
    }

    for (int j=1; j<=v.size(); j++){
      if (v(j)<0)
        v(j) = vmin;
    }
    
    return true;
}
            
bool StatModel::calcBestFitA_CG3()
{
  int i,j;
 
  // set bounds on logA
  maxlogA = -100;
  minlogA = 100;
  int kk=1;
  double diff2;

  for (i=1; i<=i_para_scaling.size(); i++){
    for (j=1; j<=data.size(); j++){
          diff2 = log(bvec(j))-log(fmat_para(j,i_para_scaling(i)));
          if (diff2>maxlogA(kk)) maxlogA(kk) = diff2;
          if (diff2<minlogA(kk)) minlogA(kk) = diff2;
    }
    kk++;
  }
  
  for (i=1; i<=i_seq_scaling.size(); i++){
    for (j=1; j<=data.size(); j++){
      diff2 = log(bvec(j))-log(fmat_seq(j,i_seq_scaling(i)));
      if (diff2>maxlogA(kk)) maxlogA(kk) = diff2;
      if (diff2<minlogA(kk)) minlogA(kk) = diff2;
    }
    kk++;
  }
  
  for (int k=1; k<=maxlogA.size(); k++){
    maxlogBbyA(k) = maxlogA(k)-minlogA(1);//subtracting min not max
    minlogBbyA(k) = minlogA(k)-maxlogA(1);//subtracting  max not min 
    if (k>1){
      maxlogBbyA(k) += 1.0;
      minlogBbyA(k) -= 1.0;
    }
    if (minlogBbyA(k) > maxlogBbyA(k)){
      double tmp = minlogBbyA(k);
      minlogBbyA(k) = maxlogBbyA(k);
      maxlogBbyA(k) = tmp;
    }
  }

#ifdef _CHECK_
  cerr << "start CG trials inside StatModel. Line 1011.\n";
#endif  
  
   double min_prev_cost = 1e99;
   for (int itrial=1; itrial<=cg_ntrial; itrial++){
      // randomly initialize logBbyA
      logBbyA(1) = 0;
      for (kk=2; kk<=logA.size(); kk++){
	 logBbyA(kk) = minlogBbyA(kk)+(maxlogBbyA(kk)-minlogBbyA(kk))*ran2(&cg_idum);
      }
      setA(logBbyA);

      double prev_cost = CG_calc_cost_and_grad3(true,grad);
      direc = -grad;

     // conjugate gradient search
      for (int iter=1; iter<=cg_iter_max; iter++){
          // minimize along given gradient
          double new_cost = CG_line_min3();

          if ((prev_cost-new_cost)<=cg_tol*data_lognorm2){
              break;
	  }
	  	  
	  // calc conjugate gradient
	  CG_calc_cost_and_grad3(true,new_grad);
	  double gg=0.0, dgg=0.0;
	  for (int i=2; i<=grad.size(); i++){
              gg += grad(i)*grad(i);
              dgg += new_grad(i)*new_grad(i);
	  }
	  if (abs(gg)<1e-10){
              break;
	  }else{
              dgg /= gg;
	  }
	  for (int i=2; i<=grad.size(); i++){
              direc(i) = -new_grad(i)+dgg*direc(i);
	  }

	  grad = new_grad;
	  prev_cost = new_cost;

      }

      for (i=2; i<=logBbyA.size(); i++){
          // check for bounds
	  if ((logBbyA(i) > maxlogBbyA(i)) || (logBbyA(i) < minlogBbyA(i))){
              prev_cost *= 10;
          }
      }
      if (prev_cost<min_prev_cost){
          min_prev_cost = prev_cost;
          bestlogBbyA = logBbyA;
      }
  }

   //logBbyA has been determined. Now,calculate A and other scaling coefficients
   for (int k=1; k<=bestlogBbyA.size(); k++){
      BbyA(k) = exp(bestlogBbyA(k)); 
   }

#ifdef _CHECK_
   cerr << "Calculated BbyA by CG StatModel.Find Bias and A.Line 1074\n";
#endif

   Array1d<double> etmp(data.size());

   if (i_para_scaling.size()>0){
     for (j=1; j<=data.size(); j++){
       etmp(j) = 0;
       int kk=1;
 
      for (i=1; i<=i_para_scaling.size(); i++){
	 etmp(j) += fmat_para(j,i_para_scaling(i))*BbyA(kk);
	 kk++;
       }
       
       if (i_seq_scaling.size()>0){
	 double val = 0;
	 for (i=1; i<=i_seq_scaling.size(); i++){
	   val += 1/(fmat_seq(j,i_seq_scaling(i))*BbyA(kk));
	   kk++;
	 }
	 etmp(j) += 1/val;
       }
     }
   }

   if (i_para_scaling.size()==0 && i_seq_scaling.size()>1){
     for (j=1; j<=data.size(); j++){
       int kk = 1;
       etmp(j) = 0;
       double val = 0.0;
       for (i=1; i<=i_seq_scaling.size(); i++){
	 val += 1/(fmat_seq(j,i_seq_scaling(i))*BbyA(kk));
	 kk++;
       }
       etmp(j) = 1/val;
     }
   }

   double a1 = 0;
   Array1d<double> tmp_Acoef(nbias);
   int k;
   int l = 0;
   int indx = 1;
   int FLAG = 1;
   double sumA = 0;
   while (FLAG==1){
     k = 0;
     a1 = 0.0;
     for (j=indx; j<=data.size(); j++){
       if (runID(j) == runID(indx)){
	 a1 += log(bvec(j))-log(etmp(j));
	 k++;
       }else break;
     }
     l++;
     indx = j;
     tmp_Acoef(l) = a1/k; // logA
     sumA += tmp_Acoef(l);
     if (j>=data.size()) FLAG=0;
   }
   
   double Acoef = sumA/l; // average logA

   if (nbias>1){
     for (l=1; l<=nbias; l++){
       bias(l)->setValue(tmp_Acoef(l) - Acoef);
     }
   }

   if (i_para_scaling.size()>0){
     int kk=1;
     for (i=1; i<=i_para_scaling.size(); i++){
       bestlogA(kk) = bestlogBbyA(kk) + Acoef;
       kk++;
     }
     
     if (i_seq_scaling.size()>0){
       for (i=1; i<=i_seq_scaling.size(); i++){
	 bestlogA(kk) = bestlogBbyA(kk) + Acoef;
	 kk++;
       }
     }
   }
    
   if (i_para_scaling.size()==0 && i_seq_scaling.size()>1){
     int kk=1;
     for (i=1; i<=i_seq_scaling.size(); i++){
       bestlogA(kk) = bestlogBbyA(kk) + Acoef;
       kk++;
     }
   }

   setA(bestlogA);
   return true;
}

bool StatModel::calcBestFitA_CG2()
{
  int i,j;
 
  // set bounds on logA
  maxlogA = -100;
  minlogA = 100;
  int kk=1;
  double diff2;

  for (i=1; i<=i_para_scaling.size(); i++){
    for (j=1; j<=data.size(); j++){
          diff2 = log(bvec(j))-log(fmat_para(j,i_para_scaling(i)));
          if (diff2>maxlogA(kk)) maxlogA(kk) = diff2;
          if (diff2<minlogA(kk)) minlogA(kk) = diff2;
    }
    kk++;
  }
  
  for (i=1; i<=i_seq_scaling.size(); i++){
    for (j=1; j<=data.size(); j++){
      diff2 = log(bvec(j))-log(fmat_seq(j,i_seq_scaling(i)));
      if (diff2>maxlogA(kk)) maxlogA(kk) = diff2;
      if (diff2<minlogA(kk)) minlogA(kk) = diff2;
    }
    kk++;
  }

  for (int k=1; k<=maxlogA.size(); k++){
      maxlogA(k) += 1.0;
      minlogA(k) -= 1.0;
  }

   double min_prev_cost = 1e99;
   for (int itrial=1; itrial<=cg_ntrial; itrial++){
      // randomly initialize logA
      for (kk=1; kk<=logA.size(); kk++){
	 logA(kk) = minlogA(kk)+(maxlogA(kk)-minlogA(kk))*ran2(&cg_idum);
      }
      setA(logA);

      double prev_cost = CG_calc_cost_and_grad3(true,grad);
      direc = -grad;

     // conjugate gradient search
      for (int iter=1; iter<=cg_iter_max; iter++){
          // minimize along given gradient
          double new_cost = CG_line_min2();

          if ((prev_cost-new_cost)<=cg_tol*data_lognorm2){
              break;
	  }
	  	  
	  // calc conjugate gradient
	  CG_calc_cost_and_grad3(true,new_grad);
	  double gg=0.0, dgg=0.0;
	  for (int i=1; i<=grad.size(); i++){
              gg += grad(i)*grad(i);
              dgg += new_grad(i)*new_grad(i);
	  }
	  if (abs(gg)<1e-10){
              break;
	  }else{
              dgg /= gg;
	  }
	  for (int i=1; i<=grad.size(); i++){
              direc(i) = -new_grad(i)+dgg*direc(i);
	  }

	  grad = new_grad;
	  prev_cost = new_cost;

      }

      for (i=1; i<=logA.size(); i++){
          // check for bounds
	  if ((logA(i) > maxlogA(i)) || (logA(i) < minlogA(i))){
              prev_cost *= 10;
          }
      }
      if (prev_cost<min_prev_cost){
          min_prev_cost = prev_cost;
          bestlogA = logA;
      }
  }

   Array1d<double> etmp(data.size());

   for (j=1; j<=data.size(); j++){
     etmp(j) = 0;

     if (i_para_noscaling.size()>0){
       for (i=1; i<=i_para_noscaling.size(); i++){
	 int ii = i_para_noscaling(i);
	 double Afix = para(ii)->scaling();
	 etmp(j) += fmat_para(j,ii)*Afix;
       }
     }

     if (i_para_scaling.size()>0){
       int kk=1;
       for (i=1; i<=i_para_scaling.size(); i++){
	 etmp(j) += fmat_para(j,i_para_scaling(i))*exp(bestlogA(kk));
	 kk++;
       }

       if (i_seq_noscaling.size()>0 || i_seq_scaling.size()>0){
	 double val = 0;
	 if (i_seq_noscaling.size()>0){
	   for (i=1; i<=i_seq_noscaling.size(); i++){
	     int ii = i_seq_noscaling(i);
	     double Afix = seq(ii)->scaling();
	     val += 1/(fmat_seq(j,ii)*Afix);
	   }
	 }
     
	 if (i_seq_scaling.size()>0){
	   for (i=1; i<=i_seq_scaling.size(); i++){
	     val += 1/(fmat_seq(j,i_seq_scaling(i))*exp(bestlogA(kk)));
	     kk++;
	   }
	 }
	 etmp(j) += 1/val;
       }
     }else{
       double val = 0;
       if (i_seq_noscaling.size()>0){
	 for (i=1; i<=i_seq_noscaling.size(); i++){
	   int ii = i_seq_noscaling(i);
	   double Afix = seq(ii)->scaling();
	   val += 1/(fmat_seq(j,ii)*Afix);
	 }
       }
       int kk=1;
       for (i=1; i<=i_seq_scaling.size(); i++){
	 val += 1/(fmat_seq(j,i_seq_scaling(i))*exp(bestlogA(kk)));
	 kk++;
       }
       etmp(j) += 1/val;
     }
   }
 
   if (nbias>1){
     double x = 0;
     int k, l = 0;
     int indx = 1;
     int FLAG = 1;
     double sumA = 0;
     while (FLAG==1){
       k = 0;
       x = 0.0;
       for (j=indx; j<=data.size(); j++){
	 if (runID(j) == runID(indx)){
	   x += log(bvec(j))-log(etmp(j));
	   k++;
	 }else break;
       }
       l++;
       bias(l)->setValue(x/k);
       indx = j;
    
       if (j>=data.size()) FLAG=0;
     }  
   }

   setA(bestlogA);
   return true;
}

void StatModel::setA(const Array1d<double>& logA_)
{
    int kk=1;
    for (int i=1; i<=i_para_scaling.size(); i++){
	para(i_para_scaling(i))->setScaling(exp(logA_(kk)));
	kk++;
    }
    for (int i=1; i<=i_seq_scaling.size(); i++){
	seq(i_seq_scaling(i))->setScaling(exp(logA_(kk)));
	kk++;
    }
}

double StatModel::CG_calc_cost_and_grad3(bool needGrad,
					Array1d<double>& grad)
{
    double cost=0.0;
    if (needGrad) grad=0.0;
    bool do_pair;

#ifdef _CHECK_    
    int k=0;
#endif

    for (int j=1; j<=data.size()-1; j++){
      int jj = j+1;
      if (jj>data.size()) break;
      while (runID(j)==runID(jj)){

         if (!less_pairs){ 
	   do_pair = true;
         }else{
             do_pair = false;
	     if (newID(j)==newID(jj)){
                do_pair = true;
             }else{
                if (j==1 || newID(j-1)!=newID(j)){
                  if (newID(jj-1)!=newID(jj)) do_pair=true;
		  if (jj==data.size() || newID(jj+1)!=newID(jj)) do_pair=true;
                }
		if (newID(j+1)!=newID(j)){
                  if (newID(jj-1)!=newID(jj)) do_pair=true;
                  if (jj==data.size() || newID(jj+1)!=newID(jj)) do_pair=true;
                }
             }
           }

	 if (do_pair){

          double edot_total1=0.0;
	  double edot_total2=0.0;
	  double edot_total_seq1=0.0;
	  double edot_total_seq2=0.0;
	  double edot_total_seq_sq1, edot_total_seq_sq2;

#ifdef _CHECK_
	  k++; // to count the number of pairs being evaluated.
#endif
	  if (para.size()>0){
	    if (i_para_noscaling.size()>0){
	      for (int i=1; i<=i_para_noscaling.size(); i++){
		int ii = i_para_noscaling(i);
		double Afix = para(ii)->scaling();
		edot_total1 += fmat_para(j,ii)*Afix;
		edot_total2 += fmat_para(jj,ii)*Afix;
	      }
	    }
	    if (i_para_scaling.size()>0){
	      for (int i=1; i<=i_para_scaling.size(); i++){
		int ii = i_para_scaling(i); 
		edot_total1 += para(ii)->predict(fmat_para(j,ii));
		edot_total2 += para(ii)->predict(fmat_para(jj,ii));
	      }
	    }
	  }
	  if (seq.size()>0){
	    double tmp1=0.0;
	    double tmp2=0.0;
	    if (i_seq_noscaling.size()>0){
	      for (int i=1; i<=i_seq_noscaling.size(); i++){
		int ii = i_seq_noscaling(i);
		double Afix = seq(ii)->scaling();
		tmp1 += 1.0/(fmat_seq(j,ii)*Afix);
		tmp2 += 1.0/(fmat_seq(jj,ii)*Afix);
	      }
	    }
	    if (i_seq_scaling.size()>0){
	      for (int i=1; i<=i_seq_scaling.size(); i++){
		int ii = i_seq_scaling(i);
		tmp1 += 1.0/seq(ii)->predict(fmat_seq(j,ii));
		tmp2 += 1.0/seq(ii)->predict(fmat_seq(jj,ii));
	      }
	    }
	    
	    edot_total_seq1 = 1.0/tmp1;
	    edot_total_seq2 = 1.0/tmp2;
	    edot_total_seq_sq1 = edot_total_seq1*edot_total_seq1;
	    edot_total_seq_sq2 = edot_total_seq2*edot_total_seq2;
	    
	    edot_total1 += edot_total_seq1;
	    edot_total2 += edot_total_seq2;
	  }

	  double misfit = 0.0;
	  
	  double num1, num2, num, den1, den2, den;
	  num = log(data(jj).strainRate()/data(j).strainRate());
	  num -= log(edot_total2/edot_total1);

	  den1 = data(j).strainRateRelError()*data(j).strainRateRelError();
	  den2 = data(jj).strainRateRelError()*data(jj).strainRateRelError();
	  den = den1 + den2;
	  misfit = num;
	  cost += (misfit*misfit)/den;

	  // calculate gradients
	  /* 1. if (i_para_scaling>0) --> (a) A*f1 + B*g1 + C*h1 +....
	                                 (b) A*f1 + 1/(1/A_f1_ + 1/B_g1_ +..)
	     2. if (i_seq_scaling>0 && i_para_scaling==0) -->
	                            (c) 1/(1/A_f1_ + 1/B_g1_ + 1/C_h1_ +..) */
	  // we don't need to calculate grad(kk) for the 1st flow law for cases (a) & (b)
	  // we don't need to calculate it for case (c) too if i_seq_noscaling.size()==0
	  // otherwise, we need to calculate grad(kk) for the 1st flow law 
	  if (needGrad){
	    if (para.size()>0){
	      int kk, i1;
	      if (i_para_noscaling.size()>0 || i_seq_noscaling.size()>0){
		kk=2; i1=2;
	      }else{
		kk=1; i1=1;
	      }
	      if (i_para_scaling.size()>1){
		for (int i=i1; i<=i_para_scaling.size(); i++){
		  double v1 = fmat_para(j,i_para_scaling(i));
		  double v2 = fmat_para(jj,i_para_scaling(i));
		  grad(kk) += misfit * ((v2/edot_total2) - (v1/edot_total1));
		  kk++;
		}
	      }
	      
	      if (i_seq_scaling.size()>0){
		for (int i=1; i<=i_seq_scaling.size(); i++){
		  int ii = i_seq_scaling(i);
		  double tmp1 = 1.0/seq(ii)->predict(fmat_seq(j,ii));
		  double tmp2 = 1.0/seq(ii)->predict(fmat_seq(jj,ii));
		  double v3 = fmat_seq(j,ii) *edot_total_seq_sq1*tmp1*tmp1;
		  double v4 = fmat_seq(jj,ii) *edot_total_seq_sq2*tmp2*tmp2;
		  grad(kk) += misfit * ((v4/edot_total2) - (v3/edot_total1));
		  kk++;
		}
	      }
	      
	    }else{ 
	      // case where para.size()==0 : only seq. flow laws 
	      if (i_seq_scaling.size()>1){
		int kk, i1;
		if (i_seq_noscaling.size()==0){
		  kk=2; i1=2;
		}else{
		  kk=1; i1=1;
		}
		for (int i=i1; i<=i_seq_scaling.size(); i++){
		  int ii = i_seq_scaling(i);
		  double tmp1 = 1.0/seq(ii)->predict(fmat_seq(j,ii));
		  double tmp2 = 1.0/seq(ii)->predict(fmat_seq(jj,ii));
		  //double v3 = fmat_seq(j,ii) *edot_total_seq_sq1*tmp1*tmp1;
		  //double v4 = fmat_seq(jj,ii) *edot_total_seq_sq2*tmp2*tmp2;
		  //grad(kk) += misfit * ((v4/edot_total2) - (v3/edot_total1));
		  // simplify this-
		  double v3 = fmat_seq(j,ii)*tmp1*tmp1;
		  double v4 = fmat_seq(jj,ii)*tmp2*tmp2;
		  grad(kk) += misfit * ((v4*edot_total2) - (v3*edot_total1));
		  kk++;
		}
	      }
	    }
	  }
	 }
	 jj++;
	  
	 if (jj>data.size()) break;
	 if (runID(j)!=runID(jj)) break;
      }
    }
    
    // finish up gradients
    if (needGrad){
      if (para.size()>0){
	int kk, i1;
	if (i_para_noscaling.size()==0 || i_seq_noscaling.size()>0){
	  kk=2; i1=2;
	}else{
	  kk=1; i1=1;
	}
	if (i_para_scaling.size()>1){
	  for (int i=i1; i<=i_para_scaling.size(); i++){
	    grad(kk) *= (-2.0)*para(i_para_scaling(i))->scaling();
	    kk++;
	  }
	}
	for (int i=1; i<=i_seq_scaling.size(); i++){
	  grad(kk) *= (-2.0)*seq(i_seq_scaling(i))->scaling();
	  kk++;
	}
      }else{
	if (i_seq_scaling.size()>1){
	  int kk, i1;
	  if (i_seq_noscaling.size()==0){
	    kk=2; i1=2;
	  }else{
	    kk=1; i1=1;
	  }
	  for (int i=i1; i<=i_seq_scaling.size(); i++){
	    grad(kk) *= (-2.0)*seq(i_seq_scaling(i))->scaling();
	    kk++;
	  }
	}
      }
    }

#ifdef _CHECK_
      cerr << "Number of pairs evaluated during CG: " << k << "\n";
#endif

    return cost;
}

double StatModel::CG_line_min3()
{
    double ax = 0.0;
    double xx = 1.0;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &StatModel::f1dim3); 
    double new_val = brent(ax,xx,bx,&xmin, &StatModel::f1dim3);
    for (int i=2; i<=logBbyA.size(); i++){
	logBbyA(i) += xmin*direc(i);
    }
    setA(logBbyA);
    return new_val;
}

double StatModel::CG_line_min2()
{
    double ax = 0.0;
    double xx = 1.0;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &StatModel::f1dim2); 
    double new_val = brent(ax,xx,bx,&xmin, &StatModel::f1dim2);
    for (int i=1; i<=logA.size(); i++){
	logA(i) += xmin*direc(i);
    }
    setA(logA);
    return new_val;
}

double StatModel::f1dim3(double x)
{
    for (int i=1; i<=logBbyA.size(); i++){
      if (i==1) tmp_logBbyA(1) = 0;
      else  tmp_logBbyA(i) = logBbyA(i)+x*direc(i);
    }
    setA(tmp_logBbyA);
    double tmp = CG_calc_cost_and_grad3(false,grad);
    return tmp;
}

double StatModel::f1dim2(double x)
{
    for (int i=1; i<=logA.size(); i++){
      tmp_logA(i) = logA(i)+x*direc(i);
    }
    setA(tmp_logA);
    double tmp = CG_calc_cost_and_grad3(false,grad);
    return tmp;
}

void StatModel::afterBestFitA3()
{

  //cerr << "dif  " << "gbs\n";
 for (int j=1; j<=data.size(); j++){
    if (i_para_noscaling.size()>0){
      for (int i=1; i<=i_para_noscaling.size(); i++){
	int ii = i_para_noscaling(i);
	//    edot_para(j,ii) = fmat_para(j,ii);
	edot_para_tmp(j,ii) = fmat_para(j,ii)*para(ii)->scaling();
	//cerr << edot_para_tmp(j,ii) << " ";
      }
    }

    int kk=1;
    if (i_para_scaling.size()>0){
      for (int i=1; i<=i_para_scaling.size(); i++){
	int ii = i_para_scaling(i);
      //edot_para(j,ii) 
      //= para(ii)->predict(fmat_para(j,ii));
	//edot_para_tmp(j,ii) = BbyA(kk) * fmat_para(j,ii);
	edot_para_tmp(j,ii) = para(ii)->predict(fmat_para(j,ii));
	//cerr << edot_para_tmp(j,ii) << "\n";
	kk++;
      }
    }

    if (i_seq_noscaling.size()>0){
      for (int i=1; i<=i_seq_noscaling.size(); i++){
	int ii = i_seq_noscaling(i);
	//  edot_seq(j,ii) = fmat_seq(j,ii);
	edot_seq_tmp(j,ii) = fmat_seq(j,ii)*seq(ii)->scaling();
      }
    }

    if (i_seq_scaling.size()>0){
      for (int i=1; i<=i_seq_scaling.size(); i++){
	int ii = i_seq_scaling(i);
	//edot_seq(j,ii) = seq(ii)->predict(fmat_seq(j,ii));
	//edot_seq_tmp(j,ii) = BbyA(kk) * fmat_seq(j,ii);
	edot_seq_tmp(j,ii) = seq(ii)->predict(fmat_seq(j,ii));
	kk++;
      }
    }
  }
}

/*void StatModel::afterBestBbyA()
{
  if (i_para_noscaling.size()>0){
    for (int i=1; i<=i_para_noscaling.size(); i++){
	int ii = i_para_noscaling(i);
	for (int j=1; j<=data.size(); j++){
    	    edot_para_tmp(j,ii) = fmat_para(j,ii);
	}
    }
  }

  int kk=1;
  if (i_para_scaling.size()>0){
    for (int i=1; i<=i_para_scaling.size(); i++){
	int ii = i_para_scaling(i);
	for (int j=1; j<=data.size(); j++){
	  edot_para_tmp(j,ii) = BbyA(kk) * fmat_para(j,ii);
	}
	kk++;
    }
  }

  if (i_seq_noscaling.size()>0){
    for (int i=1; i<=i_seq_noscaling.size(); i++){
	int ii = i_seq_noscaling(i);
	for (int j=1; j<=data.size(); j++){
	  edot_seq_tmp(j,ii) = fmat_seq(j,ii);
	}
    }
  }

  if (i_seq_scaling.size()>0){
    for (int i=1; i<=i_seq_scaling.size(); i++){
	int ii = i_seq_scaling(i);
	for (int j=1; j<=data.size(); j++){
	  edot_seq_tmp(j,ii) = BbyA(kk) * fmat_seq(j,ii);
	}
	kk++;
    }
  }
  }*/

double StatModel::calcChiSq3(int k, double& chi2_orig)
{
    fixParams();

    if (k<0 || !unfixed_params(k)->isScaling()){
        beforeBestFitA(k);
        if (calcBestFitA3() == false){

	  // comment on Jan 8 cerr << "calcFitA3 is false, chi2=" << max_chi2 << "\n";

            return max_chi2;
        }
    }
    afterBestFitA3();
    
    double chi2=0.0;
    chi2_orig=0.0;
    for (int j=1; j<=data.size()-1; j++){
        int jj = j+1;
        while (runID(j)==runID(jj)){
	 
	    double edot_total1=0.0;
            double edot_total_seq1=0.0;
            double edot_total_seq_sq1=0.0;
            double edot_total2=0.0;
            double edot_total_seq2=0.0;
            double edot_total_seq_sq2=0.0;
	    
	    if (para.size()>0){
	      for (int i=1; i<=para.size(); i++){
                edot_total1 += edot_para_tmp(j,i);
                edot_total2 += edot_para_tmp(jj,i);
	      }
            }

            if (seq.size()>0){
                double tmp1=0.0;
                double tmp2=0.0;
                for (int i=1; i<=seq.size(); i++){
                    tmp1 += 1.0/edot_seq_tmp(j,i);
                    tmp2 += 1.0/edot_seq_tmp(jj,i);
                }
                edot_total_seq1 = 1.0/tmp1;
                edot_total_seq2 = 1.0/tmp2;
                edot_total_seq_sq1 = edot_total_seq1*edot_total_seq1;
                edot_total_seq_sq2 = edot_total_seq2*edot_total_seq2;
                edot_total1 += edot_total_seq1;
                edot_total2 += edot_total_seq2;
            }
            
            double den = 0.0;
            double v11 = data(j).strainRateRelError();
            double v12 = data(jj).strainRateRelError();
            double val0sq = v11*v11 + v12*v12;
            den += val0sq;
	    
            if (para.size()>0){
	      for (int i=1; i<=para.size(); i++){
                double v21 = para(i)->variance(data(j));
                double v22 = para(i)->variance(data(jj));
                double v31 = edot_para_tmp(j,i)/edot_total1;
                double v32 = edot_para_tmp(jj,i)/edot_total2;
                double v41 = v31*v31*v21;
                double v42 = v32*v32*v22;
                den += v41 + v42;
	      }
	    }

	    if (seq.size()>0){
	      for (int i=1; i<=seq.size(); i++){
                double v51 = edot_total_seq_sq1/(edot_total1*edot_seq_tmp(j,i));
                double v52 = edot_total_seq_sq2/(edot_total2*edot_seq_tmp(jj,i));
                double v61 = v51*v51*seq(i)->variance(data(j));
                double v62 = v52*v52*seq(i)->variance(data(jj));
                den += v61 + v62;
	      }
	    }

            double misfit = log(data(jj).strainRate()) - log(data(j).strainRate());
	    misfit -= (log(edot_total2) - log(edot_total1));
            double num = misfit*misfit;

	    //            chivec(j) = num/val0sq;
            chi2 += num/val0sq; // simple chi2
            chi2_orig += num/den; // full chi2
	    jj++;
	    if (jj>data.size()) break;
	    if (runID(j)!=runID(jj)) break;
	}
    }

    return chi2;
}

