//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "aux_histogram.hxx"
#include "aux_utils.hxx"
#include "mprintf.hxx"
#include "Mesh/MeshBase.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include <unistd.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <cstdio>
#include <algorithm>

namespace Metris{


void print_histogram(const MeshBase &msh, dblAr1 &values, IntrpTyp iinter,
                     dblAr1 &bounds, std::string symb, std::string name){
  GETVDEPTH(msh);
  if(!DOPRINTS1()) return;

  int nval = values.get_n();
  if(nval <= 0) return;


  const int nbucket = 10; 

  dblAr2 buckval(nbucket,2);
  intAr1 buckcnt(nbucket);
  buckval.set_n(nbucket);
  buckcnt.set_n(nbucket);
  buckcnt.fill(nbucket, 0);
  double vlow, vhig;

  double vmin = 1.0e30, vmax = -1.0e30, vavgl = 0, vavgg = 0;

  if(bounds.get_n() == 0){
    vlow = values[0];
    vhig = values[0];
    std::sort(&values[0], &values[nval], std::less<double>{});
    // The 1/nbucket smallest and largest are excluded -> often outliers. 
    vlow = values[nval / nbucket];
    vhig = values[nval - nval / nbucket];
  }else{
    METRIS_ASSERT(bounds.get_n() == 2);
    vlow = bounds[0];
    vhig = bounds[1];
  }



  buckval(0,0) = vlow;
  for(int ibucket = 0; ibucket < nbucket-1; ibucket++){
    double val; 
    if(iinter == IntrpTyp::Linear){
      val = vlow + (ibucket+1) * (vhig - vlow) / nbucket;
    }else{
      val = vlow * pow(vhig/vlow,(double)(ibucket+1) / (double)nbucket);
    }
    buckval(ibucket,1) = val;
    if(ibucket > 0) buckval(ibucket,0) = buckval(ibucket-1,1);
  }

  buckval(nbucket-1,0) = buckval(nbucket-2,1);
  buckval(nbucket-1,1) = vhig;

  int nlow = 0;
  int nhig = 0;
  bool nogeom = false;
  int imax = -1;
  int imin = -1;
  int ival = 0;
  for(double val : values){
    if(val < vmin){
      vmin = val;
      imin = ival;
    }
    if(val > vmax){
      vmax = val;
      imax = ival;
    }
    vavgl += val;
    if(val < 1.0e-16) nogeom = true;
    if(!nogeom) vavgg += log(val);

    for(int ibucket = 0; ibucket < nbucket; ibucket++){
      if(val < buckval(ibucket,1) && val >= buckval(ibucket,0)){
        buckcnt[ibucket]++;
      }
    }
    if(val < vlow) nlow++;
    if(val > vhig) nhig++;
    ival++;
  }

  if(!nogeom) vavgg = exp(vavgg / nval);
  vavgl = vavgl / nval;

  int maxbuckt = 0;
  for(int nn : buckcnt) maxbuckt = MAX(maxbuckt, nn);
  maxbuckt = MAX(nlow,maxbuckt);
  maxbuckt = MAX(nhig,maxbuckt);


  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  const int mcol = 100;
  w.ws_col = MIN(MAX(w.ws_col, 40),mcol);

  // Expected length to print %8.2e < symb < %8.2e w/ spaces
  int ibuf0  = symb.length() + 6 + 2*8 + 2; 
  int prtwdt = w.ws_col - ibuf0;

  // The maximum should not exceed height
  // This determines the scaling factor 
  double scal = prtwdt / (double) maxbuckt; 
  
  int ncol = MIN(mcol,MAX(w.ws_col, ibuf0 + 15));
  char buffer[2+nbucket][mcol];

  for(int ii = 0; ii < nbucket; ii++){
    for(int jj = 0; jj < ncol; jj++){
      buffer[1+ii][jj] = ' ';
    }
  }

  CPRINTF1("-- %s: %.2f %% within bounds %f %f \n", name.c_str(),
                           100.0 - 100.0*(nlow + nhig)/(double) nval,vlow, vhig);
  CPRINTF1("  - minimum = %f (%d)\n",vmin,imin);
  CPRINTF1("  - maximum = %f (%d)\n",vmax,imax);
  if(!nogeom){
    CPRINTF1("  - average = %f (geometric) = %f\n",vavgl,vavgg);
  }else{
    CPRINTF1("  - average = %f \n",vavgl);
  }

  if(!DOPRINTS2()) return;

  // Bucket labels: 
  if(iinter == IntrpTyp::Linear){
    snprintf(buffer[0],ibuf0,"           %s < %8.2f :",symb.c_str(),buckval(0,0));
  }else{
    snprintf(buffer[0],ibuf0,"           %s < %8.2e :",symb.c_str(),buckval(0,0));
  }
  for(int ii = 0; ii < nbucket; ii++){
    if(iinter == IntrpTyp::Linear){
      snprintf(buffer[1+ii],ibuf0,"%8.2f < %s < %8.2f :",
                                      buckval(ii,0),symb.c_str(),buckval(ii,1));
    }else{
      snprintf(buffer[1+ii],ibuf0,"%8.2e < %s < %8.2e :",
                                      buckval(ii,0),symb.c_str(),buckval(ii,1));
    }
  }
  if(iinter == IntrpTyp::Linear){
    snprintf(buffer[nbucket+1],ibuf0,"%8.2f < %s            :",buckval(nbucket-1,1),
                                                         symb.c_str());
  }else{
    snprintf(buffer[nbucket+1],ibuf0,"%8.2e < %s            :",buckval(nbucket-1,1),
                                                         symb.c_str());
  }

  // Bucket lines: 
  for(int ibucket = -1; ibucket <= nbucket; ibucket++){
    int nchar;
    if(ibucket == -1){
      nchar = (int) (scal * nlow); 
    }else if(ibucket == nbucket){
      nchar = (int) (scal * nhig); 
    }else{
      nchar = (int) (scal * buckcnt[ibucket]);
    }
    METRIS_ASSERT_MSG(nchar <= ncol,"nchar = "<<nchar<<" ncol = "<<ncol);
    for(int ii = 0; ii < nchar; ii++){
      buffer[1+ibucket][ibuf0 + ii] = '*';
    }
    for(int ii = nchar; ii < ncol - ibuf0; ii++){
      buffer[1+ibucket][ibuf0 + ii] = ' ';
    }
  }
 
  for(int ii = 0; ii < nbucket+2; ii++){
    for(int jj = 0; jj < ncol; jj++){
      printf("%c",buffer[ii][jj]);
    }
    printf("\n");
  }

}



}//end namespace