/*
 * flowlaw_new.cc
 *
 * Jun Korenaga
 * Summer 2008
 *
 * CJ: 
 * Last modified in June'17, Sept.'17
 * added a fourth flow law for Peierls mechanism (Peierls4)
 * added a gbs wet case
 * Februry 2018: included provision for fixed scaling coeff.
 */

#include <cmath>
#include "constants.h"
#include "flowlaw.h"
#include "labmc_nr.h"

// FlowLaw
bool FlowLaw::checkStates(const Array1d<int>& states) const
{
    int ncount=0;
    for (int i=1; i<=required_states.size(); i++){
	int ii=required_states(i);
	for (int j=1; j<=states.size(); j++){
	    if (states(j) == ii){
		ncount++;
		continue;
	    }
	}
    }

    return (ncount == required_states.size()) ? true : false;
}

//FlowLawGen
FlowLawGen::FlowLawGen(	double facA_min, double facA_max,
		   	double n_min, double n_max,
		   	double s_min, double s_max,
		   	double E_min, double E_max)
{
    required_states.resize(1);
    required_states(1) = State::oxygen_fugacity;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Gen:facA",true);
	unfixed_params.push_back(facA__);
	// Note: even though facA__ is unfixed, I don't count it
	// among unfixed_params simply because scaling coeff. are
	// treated differently from other flow-law parameters in
	// our inversion scheme.
	// in the current inversion scheme, this condition is never met.
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"Gen:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  A = pow(10,facA_min);
	  unfixedScaling = false;
	}
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"Gen:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"Gen:n");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Gen:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Gen:E");
    }else{
	nerror++;
    }

    if (s_min < s_max){
	s__ = new Parameter(s_min,s_max,"Gen:s");
	unfixed_params.push_back(s__);
    }else if (s_min == s_max){
	s__ = new Parameter(s_min,"Gen:s");
    }else{
	nerror++;
    }
    
    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawGen::~FlowLawGen()
{
    delete facA__;
    delete n__;
    delete E__;
    delete s__;
}

double FlowLawGen::prepPredict(const State& s) const
{
    double val1 = pow(s.oxygenFugacity(), s__->value());
    double val2 = pow(s.stress(), n__->value());
    double val3 = E__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*exp(-val3/val4);
}

double FlowLawGen::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
 
}

double FlowLawGen::predict(double prep) const
{
   return prep*A*exp(facA__->value());
}

double FlowLawGen::variance(const State& s) const
{
    double val1 = (s__->value())*s.oxygenFugacityRelError();
    double val2 = (n__->value())*s.stressRelError();
    double T = s.temperature();

    double val3
	= ((E__->value())/(Rgas*T*T))
	*s.temperatureError();


    double var = val1*val1+val2*val2+val3*val3;
    return var;
}


// FlowLawDiffDry
FlowLawDiffDry::FlowLawDiffDry(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::grain_size;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DiffDry:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"DiffDry:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	    unfixedScaling = false;
	    A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"DiffDry:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"DiffDry:m");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DiffDry:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DiffDry:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DiffDry:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DiffDry:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDiffDry::~FlowLawDiffDry()
{
    delete facA__;
    delete m__;
    delete E__;
    delete V__;
}

double FlowLawDiffDry::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = E__->value()+s.pressure()*V__->value();
    double val3 = Rgas*s.temperature();

    return val1*s.stress()*exp(-val2/val3);
}

double FlowLawDiffDry::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDiffDry::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDiffDry::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = s.stressRelError();
    double T = s.temperature();
    double val3 = (V__->value()/(Rgas*T))*s.pressureError();
    double val4
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4;
    return var;
}

// FlowLawDisDry
FlowLawDisDry::FlowLawDisDry(double facA_min, double facA_max,
			     double n_min, double n_max,
			     double E_min, double E_max,
			     double V_min, double V_max)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DisDry:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"DisDry:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"DisDry:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"DisDry:n");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DisDry:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DisDry:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DisDry:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DisDry:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDisDry::~FlowLawDisDry()
{
    delete facA__;
    delete n__;
    delete E__;
    delete V__;
}

double FlowLawDisDry::prepPredict(const State& s) const
{
    double val1 = pow(s.stress(), n__->value());
    double val2 = E__->value()+s.pressure()*V__->value();
    double val3 = Rgas*s.temperature();

    return val1*exp(-val2/val3);
}

double FlowLawDisDry::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDisDry::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDisDry::variance(const State& s) const
{
    double val1 = (n__->value())*s.stressRelError();
    double T = s.temperature();
    double val2 = (V__->value()/(Rgas*T))*s.pressureError();
    double val3
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3;
    return var;
}

// FlowLawDiffWet
FlowLawDiffWet::FlowLawDiffWet(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double r_min, double r_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(2);
    required_states(1) = State::grain_size;
    required_states(2) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DiffWet:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"DiffWet:facA",true);
	if (facA_min==0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"DiffWet:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"DiffWet:m");
    }else{
	nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"DiffWet:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"DiffWet:r");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DiffWet:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DiffWet:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DiffWet:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DiffWet:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDiffWet::~FlowLawDiffWet()
{
    delete facA__;
    delete m__;
    delete r__;
    delete E__;
    delete V__;
}

double FlowLawDiffWet::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = pow(s.waterContent(), r__->value());
    double val3 = E__->value()+s.pressure()*V__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*s.stress()*exp(-val3/val4);
}

double FlowLawDiffWet::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDiffWet::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDiffWet::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = (r__->value())*s.waterContentRelError();
    double val3 = s.stressRelError();
    double T = s.temperature();
    double val4 = (V__->value()/(Rgas*T))*s.pressureError();
    double val5
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4+val5*val5;
    return var;
}

// FlowLawDisWet
FlowLawDisWet::FlowLawDisWet(double facA_min, double facA_max,
			     double n_min, double n_max,
			     double r_min, double r_max,
			     double E_min, double E_max,
			     double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"DisWet:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"DisWet:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"DisWet:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"DisWet:r");
    }else{
	nerror++;
    }

    if (n_min < n_max){
	n__ = new Parameter(n_min,n_max,"DisWet:n");
	unfixed_params.push_back(n__);
    }else if (n_min == n_max){
	n__ = new Parameter(n_min,"DisWet:n");
    }else{
 	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"DisWet:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"DisWet:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"DisWet:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"DisWet:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawDisWet::~FlowLawDisWet()
{
    delete facA__;
    delete n__;
    delete E__;
    delete V__;
}

double FlowLawDisWet::prepPredict(const State& s) const
{
    double val1 = pow(s.waterContent(), r__->value());
    double val2 = pow(s.stress(), n__->value());
    double val3 = E__->value()+s.pressure()*V__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*exp(-val3/val4);
}

double FlowLawDisWet::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawDisWet::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawDisWet::variance(const State& s) const
{
    double val1 = (r__->value())*s.waterContentRelError();
    double val2 = (n__->value())*s.stressRelError();
    double T = s.temperature();
    double val3 = (V__->value()/(Rgas*T))*s.pressureError();
    double val4
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4;
    return var;
}

// FlowLawPeierls
FlowLawPeierls::FlowLawPeierls(double facA_min, double facA_max,
			       double sigP0_min, double sigP0_max,
			       double E_min, double E_max,
			       double V_min, double V_max,
			       double p, double q)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Peierls:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"Peierls:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (sigP0_min < sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,sigP0_max*1e9,"Peierls:sigP0");
	unfixed_params.push_back(sigP0__);
    }else if (sigP0_min == sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,"Peierls:sigP0");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Peierls:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Peierls:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"Peierls:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"Peierls:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
    q1 = p; q2 = q;
}

FlowLawPeierls::~FlowLawPeierls()
{
    delete facA__;
    delete sigP0__;
    delete E__;
    delete V__;
}

double FlowLawPeierls::prepPredict(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val2 = 1+(dGp*s.pressure()/G0);
  double val3 = pow((val1/val2), q1);
  double val4 = pow((1-val3), q2);
  double val5 = E__->value()+s.pressure()*V__->value();
  double val6 = Rgas*s.temperature();

  return s.stress()*s.stress()*exp(-val4*val5/val6);
}

double FlowLawPeierls::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawPeierls::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawPeierls::variance(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val2 = G0+(dGp*s.pressure());
  double val3 = pow(val1*G0/val2, q1);
  double val4 = pow((1-val3), (q2-1));
  double val5 = (1-val3)*val4;
  double val6 = Rgas*s.temperature();
  double val7 = (E__->value()+s.pressure()*V__->value())/val6;
  double val8 = V__->value()*s.pressure()/val6;
  double val9 = dGp*s.pressure()/val2;

  double val11 = s.stressRelError()*(2 + q1*q2*val7*val3*val4);
  double val12 = -1*s.pressureRelError()*(val5*val8 + q1*q2*val3*val4*val7*val9);
  double val13 = s.temperatureRelError()*val7*val5;
 
  //  cerr << val11 << " " << val12 << " " << val13 << '\n';
    double var = val11*val11+val12*val12+val13*val13;
    return var;
}


// FlowLawPeierls2 : modified Peierls (w/o sigma^2 term)
FlowLawPeierls2::FlowLawPeierls2(double facA_min, double facA_max,
			       double sigP0_min, double sigP0_max,
			       double E_min, double E_max,
			       double V_min, double V_max,
			       double p, double q)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Peierls2:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"Peierls2:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (sigP0_min < sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,sigP0_max*1e9,"Peierls2:sigP0");
	unfixed_params.push_back(sigP0__);
    }else if (sigP0_min == sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,"Peierls2:sigP0");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Peierls2:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Peierls2:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"Peierls2:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"Peierls2:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
    q1 = p; q2 = q;
}

FlowLawPeierls2::~FlowLawPeierls2()
{
    delete facA__;
    delete sigP0__;
    delete E__;
    delete V__;
}

double FlowLawPeierls2::prepPredict(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val2 = 1+(dGp*s.pressure()/G0);
  double val3 = pow((val1/val2), q1);
  double val4 = pow((1-val3), q2);
  double val5 = E__->value()+s.pressure()*V__->value();
  double val6 = Rgas*s.temperature();

  return exp(-val4*val5/val6);
}

double FlowLawPeierls2::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawPeierls2::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawPeierls2::variance(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val2 = G0+(dGp*s.pressure());
  double val3 = pow(val1*G0/val2, q1);
  double val4 = pow((1-val3), (q2-1));
  double val5 = (1-val3)*val4;
  double val6 = Rgas*s.temperature();
  double val7 = (E__->value()+s.pressure()*V__->value())/val6;
  double val8 = V__->value()*s.pressure()/val6;
  double val9 = dGp*s.pressure()/val2;

  double val11 = s.stressRelError()*(q1*q2*val7*val3*val4);
  double val12 = -1*s.pressureRelError()*(val5*val8 + q1*q2*val3*val4*val7*val9);
  double val13 = s.temperatureRelError()*val7*val5;
 
  //  cerr << val11 << " " << val12 << " " << val13 << '\n';
    double var = val11*val11+val12*val12+val13*val13;
    return var;
}


// FlowLawPeierls3: Peierls modified to exclude pressure effects bu sigma^2 term remains
FlowLawPeierls3::FlowLawPeierls3(double facA_min, double facA_max,
			       double sigP0_min, double sigP0_max,
			       double E_min, double E_max,
			       double p, double q)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Peierls3:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"Peierls3:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (sigP0_min < sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,sigP0_max*1e9,"Peierls3:sigP0");
	unfixed_params.push_back(sigP0__);
    }else if (sigP0_min == sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,"Peierls3:sigP0");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Peierls3:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Peierls3:E");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
    q1 = p; q2 = q;
}

FlowLawPeierls3::~FlowLawPeierls3()
{
    delete facA__;
    delete sigP0__;
    delete E__;
}

double FlowLawPeierls3::prepPredict(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val3 = pow(val1, q1);
  double val4 = pow((1-val3), q2);
  double val5 = E__->value();
  double val6 = Rgas*s.temperature();

  return s.stress()*s.stress()*exp(-val4*val5/val6);
}

double FlowLawPeierls3::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawPeierls3::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawPeierls3::variance(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val3 = pow(val1, q1);
  double val4 = pow((1-val3), (q2-1));
  double val5 = (1-val3)*val4;
  double val6 = Rgas*s.temperature();
  double val7 = (E__->value())/val6;

  double val11 = s.stressRelError()*(2 + q1*q2*val7*val3*val4);
  double val13 = s.temperatureRelError()*val7*val5;
 
  double var = val11*val11+val13*val13;
  return var;
}


// FlowLawPeierls4: Peierls modified to exclude pressure effects & sigma^2 term
FlowLawPeierls4::FlowLawPeierls4(double facA_min, double facA_max,
			       double sigP0_min, double sigP0_max,
			       double E_min, double E_max,
			       double p, double q)
{
    required_states.resize(0);

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"Peierls4:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"Peierls4:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (sigP0_min < sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,sigP0_max*1e9,"Peierls4:sigP0");
	unfixed_params.push_back(sigP0__);
    }else if (sigP0_min == sigP0_max){
	sigP0__ = new Parameter(sigP0_min*1e9,"Peierls4:sigP0");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"Peierls4:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"Peierls4:E");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
    q1 = p; q2 = q;
}

FlowLawPeierls4::~FlowLawPeierls4()
{
    delete facA__;
    delete sigP0__;
    delete E__;
}

double FlowLawPeierls4::prepPredict(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val3 = pow(val1, q1);
  double val4 = pow((1-val3), q2);
  double val5 = E__->value();
  double val6 = Rgas*s.temperature();

  return exp(-val4*val5/val6);
}

double FlowLawPeierls4::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawPeierls4::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawPeierls4::variance(const State& s) const
{
  double val1 = s.stress()*1e6/sigP0__->value();
  double val3 = pow(val1, q1);
  double val4 = pow((1-val3), (q2-1));
  double val5 = (1-val3)*val4;
  double val6 = Rgas*s.temperature();
  double val7 = (E__->value())/val6;

  double val11 = s.stressRelError()*(q1*q2*val7*val3*val4);
  double val13 = s.temperatureRelError()*val7*val5;
 
  double var = val11*val11+val13*val13;
  return var;
}


// FlowLawGBSDry
FlowLawGBSDry::FlowLawGBSDry(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double n_min, double n_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::grain_size;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"GBSDry:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"GBSDry:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"GBSDry:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"GBSDry:m");
    }else{
	nerror++;
    }

    if (n_min < n_max){
      n__ = new Parameter(n_min,n_max,"GBSDry:n");
      unfixed_params.push_back(n__);
    }else if (n_min == n_max){
      n__ = new Parameter(n_min,"GBSDry:n");
    }else{
      nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"GBSDry:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"GBSDry:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"GBSDry:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"GBSDry:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawGBSDry::~FlowLawGBSDry()
{
    delete facA__;
    delete m__;
    delete n__;
    delete E__;
    delete V__;
}

double FlowLawGBSDry::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = pow(s.stress(), n__->value());
    double val3 = E__->value()+s.pressure()*V__->value();
    double val4 = Rgas*s.temperature();

    return val1*val2*exp(-val3/val4);
}

double FlowLawGBSDry::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawGBSDry::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawGBSDry::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = (n__->value())*s.stressRelError();
    double T = s.temperature();
    double val3 = (V__->value()/(Rgas*T))*s.pressureError();
    double val4
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4;
    return var;
}


// FlowLawGBSWet
FlowLawGBSWet::FlowLawGBSWet(double facA_min, double facA_max,
			       double m_min, double m_max,
			       double n_min, double n_max,
			       double r_min, double r_max,
			       double E_min, double E_max,
			       double V_min, double V_max)
{
    required_states.resize(1);
    required_states(1) = State::grain_size;
    required_states(2) = State::water_content;

    int nerror = 0;
    unfixed_params.resize(0);

    if (facA_min < facA_max){
	facA__ = new Parameter(facA_min,facA_max,"GBSWet:facA",true);
	unfixed_params.push_back(facA__);
	unfixedScaling = true;
	A = 1.0;
    }else if (facA_min == facA_max){
	facA__ = new Parameter(0,"GBSWet:facA",true);
	if (facA_min == 0){
	  unfixedScaling = true;
	  A = 1.0;
	}else{
	  unfixedScaling = false;
	  A = pow(10,facA_min);
	}
    }else{
	nerror++;
    }

    if (m_min < m_max){
	m__ = new Parameter(m_min,m_max,"GBSWet:m");
	unfixed_params.push_back(m__);
    }else if (m_min == m_max){
	m__ = new Parameter(m_min,"GBSWet:m");
    }else{
	nerror++;
    }

    if (n_min < n_max){
      n__ = new Parameter(n_min,n_max,"GBSWet:n");
      unfixed_params.push_back(n__);
    }else if (n_min == n_max){
      n__ = new Parameter(n_min,"GBSWet:n");
    }else{
      nerror++;
    }

    if (r_min < r_max){
	r__ = new Parameter(r_min,r_max,"GBSWet:r");
	unfixed_params.push_back(r__);
    }else if (r_min == r_max){
	r__ = new Parameter(r_min,"GBSWet:r");
    }else{
	nerror++;
    }

    if (E_min < E_max){
	E__ = new Parameter(E_min*1e3,E_max*1e3,"GBSWet:E");
	unfixed_params.push_back(E__);
    }else if (E_min == E_max){
	E__ = new Parameter(E_min*1e3,"GBSWet:E");
    }else{
	nerror++;
    }

    if (V_min < V_max){
	V__ = new Parameter(V_min*1e-6,V_max*1e-6,"GBSWet:V");
	unfixed_params.push_back(V__);
    }else if (V_min == V_max){
	V__ = new Parameter(V_min*1e-6,"GBSWet:V");
    }else{
	nerror++;
    }

    for (int i=1; i<=unfixed_params.size(); i++){
	unfixed_params(i)->randomize();
    }
}

FlowLawGBSWet::~FlowLawGBSWet()
{
    delete facA__;
    delete m__;
    delete n__;
    delete r__;
    delete E__;
    delete V__;
}

double FlowLawGBSWet::prepPredict(const State& s) const
{
    double val1 = pow(s.grainSize(), -1.0*(m__->value()));
    double val2 = pow(s.stress(), n__->value());
    double val3 = pow(s.waterContent(), r__->value());
    double val4 = E__->value()+s.pressure()*V__->value();
    double val5 = Rgas*s.temperature();

    return val1*val2*val3*exp(-val4/val5);
}

double FlowLawGBSWet::predict(const State& s) const
{
    double val1 = prepPredict(s);

    return val1*A*exp(facA__->value());
}

double FlowLawGBSWet::predict(double prep) const
{
    return prep*A*exp(facA__->value());
}

double FlowLawGBSWet::variance(const State& s) const
{
    double val1 = (m__->value())*s.grainSizeRelError();
    double val2 = (n__->value())*s.stressRelError();
    double val3 = (r__->value())*s.waterContentRelError();
    double T = s.temperature();
    double val4 = (V__->value()/(Rgas*T))*s.pressureError();
    double val5
	= ((E__->value()+s.pressure()*V__->value())/(Rgas*T*T))
	*s.temperatureError();

    double var = val1*val1+val2*val2+val3*val3+val4*val4+val5*val5;
    return var;
}
