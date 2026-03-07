/*
 * debug.h
 *
 * Chhavi Jain
 * October 2016
 */

#ifndef _CG_DEBUG_
#define _CG_DEBUG_

using namespace std; 
#include <string>
#include <iostream>
#include<fstream>

class Debug {
public:
  //  Debug();
  //  ~Debug();

  void setEntries(int, int, int);
  void debug_file(char* in_line);
  //  void setMC(int a){imcmc=a; cerr << imcmc << '\n';};
  //  int IMCMC() {return imcmc;};

private:

  int imcmc, iflag, icp;
  
  std::string fname = "debug_VX_CC5.dat";

};
#endif /* _CG_DEBUG_ */
