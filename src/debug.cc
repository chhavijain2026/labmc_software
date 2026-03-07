/*
 * debug.cc
 *
 * Chhavi Jain
 * October 2015
 */

#include <iostream>
#include <cstdlib>
#include "debug.h"

void Debug::setEntries(int ii, int flag, int iii)
{
  imcmc = ii;
  iflag = flag;
  icp = iii;

  //  char line0[512];
  //  sprintf(line0, "%3d %3d %5.3d %5.3e %5.3e %5.3f\n",ii,flag,iii,0.0,0.0,0.0);
  //  ofstream CGdebugfile;
  //  CGdebugfile.open(fname, ios::app);
  //  CGdebugfile << imcmc << " " << iflag << " " << icp;
  //  CGdebugfile << " " << 0.0 << " " << 0.0 << " " << 0.0 << '\n';
  //  CGdebugfile << line0;
  //  CGdebugfile.close();
}

void Debug::debug_file(char* in_line)
{
  ofstream CGdebugfile;
  CGdebugfile.open (fname, ios::app);
  CGdebugfile << imcmc << " " << iflag << " " << icp << " "; 
  CGdebugfile << in_line;
  CGdebugfile.close();

}

