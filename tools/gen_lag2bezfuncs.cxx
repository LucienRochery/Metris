//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../src/ho_constants.hxx"
#include "../src/low_eval.hxx"
#include "../src/aux_utils.hxx"
#include "../src/CT_loop.hxx"


#include <boost/hana.hpp> 
namespace hana = boost::hana;
using namespace hana::literals;

#include <sstream>
#include <fstream>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>

const int nszfld = 4;
const int szflds[nszfld] = {1, 2, 3, 6};


using namespace Metris;


/*
This is as hacky as it gets, and will undoubtedly fail for high METRIS_MAX_DEG_JACOBIAN
When that time comes, use proper rational approximation methods
Works up until deg 4 included
*/
int dbl2frac(double tol, double x,  long long *p,  long long *q, bool iverb = false){
  for(int i = 2; i < 1000000; i++){
    int ix = x*i > 0 ? (int)(x*i + 0.5) : (int)(x*i - 0.5);
    if(iverb && i < 100) printf("Debug ix = %d i*x = %f err = %f\n",ix,i*x,abs(ix - x*i));
    if(abs(ix - x*i) < tol*i){
      *q = i;
      *p = ix;
      //printf("Found %f = %d / %d (= %f) up to %f \n",
      //  x,*p,*q,((double)*p)/(*q),abs(ix - i*x));
      return 0;
    }
  }

  int xp = log10(abs(x));
  *p = (long long int) (x*1e15/pow(10,xp));
  *q = (long long int) (1e15/pow(10,xp));

  //printf("Last resort approximation by %lld/%lld = %30.15f vs %30.15f\n",
  //  *p,*q,((double)*p)/ *q,x);
//  printf("Debug x 560 %30.23f\n",x*560);
  return 0;
}


template<int tdim>
void gen_lag2bez(std::ostringstream &str);

template<int tdim>
void gen_lag2bez(std::ostringstream &str){
  printf("-- Gen lag2bez tdim = %d maxdeg = %d \n",tdim,METRIS_MAX_DEG_JACOBIAN);
  double bary[tdim+1];
  char i_s[16]; 
  char j_s[16]; 
  //char v_s[64];
  //char p_s[64];
  //char q_s[64];
  char ideg_s[64];

  CT_FOR0_INC(0,METRIS_MAX_DEG_JACOBIAN,ideg){

    if constexpr(ideg == 0){
      str << "template<> void lag2bez"<<tdim<<"<"<<0<<","<<1<<">"<<
                               "(const int* __restrict__ lfld,\n\
                                const dblAr2& __restrict__ rfld0,\n\
                                dblAr2& __restrict__ rfld1){}\n\n";
      return c_ideg+1_c; // equivalent of a continue statement
    }


    snprintf(ideg_s,4,"%3d",ideg);
    int npp;
    if(tdim == 1){
    	npp = edgnpps[ideg];
    }else if(tdim == 2){
    	npp = facnpps[ideg];
    }else{
    	npp = tetnpps[ideg];
    }
    
    printf("-- Lag2bez deg = %d \n",ideg);

    Eigen::MatrixXd M(npp,npp);
    //Eigen::SparseMatrix<double> M(npp,npp);
    for(int irnk = 0; irnk < npp; irnk++){
      // The Lagrange node
      if(tdim == 1){
      	bary[0] = ordedg.s[ideg][irnk][0]/(1.0*ideg);
      	bary[1] = ordedg.s[ideg][irnk][1]/(1.0*ideg);
      	// Evaluate each Bernstein at the Lagrange
      	for(int j = 0; j < npp; j++){
          double tmp;
          //M.insert(irnk,j) = tmp = eval_bezierfunc<ideg,tdim>(ordedg.s[ideg][j],bary,0,NULL);
          M(irnk,j) = tmp = eval_bezierfunc<ideg,tdim>(ordedg.s[ideg][j],bary,0,NULL);
          //printf("Debug irnk = %d bary = %f %f val = %f\n",irnk,bary[0],bary[1],tmp);
      	}
      }else if(tdim == 2){
      	bary[0] = ordfac.s[ideg][irnk][0]/(1.0*ideg);
      	bary[1] = ordfac.s[ideg][irnk][1]/(1.0*ideg);
      	bary[2] = ordfac.s[ideg][irnk][2]/(1.0*ideg);
      	// Evaluate each Bernstein at the Lagrange
      	for(int j = 0; j < npp; j++){
          //double tmp = eval_bezierfunc<ideg,tdim>(ordfac.s[ideg][j],bary,0,NULL);
          //if(j == npp-1 && ideg == 3){
          //  printf("Debug gen lag2bez funs ideg = 3, irnk = face bary = %f %f %f , eval at node %d (%d%d%d)= %f \n",
          //            bary[0],bary[1],bary[2],j,ordfac.s[ideg][j][0],ordfac.s[ideg][j][1],ordfac.s[ideg][j][2],tmp);
          //}
          //M.insert(irnk,j) = eval_bezierfunc<ideg,tdim>(ordfac.s[ideg][j],bary,0,NULL);
      	  M(irnk,j) = eval_bezierfunc<ideg,tdim>(ordfac.s[ideg][j],bary,0,NULL);
      	}
      }else{
      	bary[0] = ordtet.s[ideg][irnk][0]/(1.0*ideg);
      	bary[1] = ordtet.s[ideg][irnk][1]/(1.0*ideg);
      	bary[2] = ordtet.s[ideg][irnk][2]/(1.0*ideg);
      	bary[3] = ordtet.s[ideg][irnk][3]/(1.0*ideg);
      	// Evaluate each Bernstein at the Lagrange
      	for(int j = 0; j < npp; j++){
          //M.insert(irnk,j) = eval_bezierfunc<ideg,tdim>(ordtet.s[ideg][j],bary,0,NULL);
      	  M(irnk,j) = eval_bezierfunc<ideg,tdim>(ordtet.s[ideg][j],bary,0,NULL);
      	}
      }
    }

    //std::cout<<"Matrix M = \n"<<M<<"\n";
    
    #if 0 
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.analyzePattern(M); 
    solver.factorize(M);
    solver.solve(M);
    if(solver.info() != Eigen::Success){
      printf("Eigen vailed to factorize matrix M !\n");
      std::cout<<"info = "<<solver.info()<<"\n";
    }
    Eigen::SparseMatrix<double> Id(npp,npp);
    Id.setIdentity();
    Eigen::SparseMatrix<double> I = solver.solve(Id);
    #else
    Eigen::MatrixXd invM = M.inverse();
    //Id(npp,npp);
    //Id.setIdentity();
    //Egen::VectorXd x = M.colPivHouseholderQr().solve(Id);
    #endif 


    std::vector<std::vector<long long int>> pcoefs(npp);
    std::vector<std::vector<long long int>> qcoefs(npp);
    std::vector<std::vector<int>> idxs(npp);

    for(int i =0; i < npp; i++){
    	pcoefs[i].reserve(npp);
    	pcoefs[i].resize(0);
      qcoefs[i].reserve(npp);
      qcoefs[i].resize(0);
    	idxs[i].reserve(npp);
    	idxs[i].resize(0);
    }

    #if 0
    for (int k=0; k<I.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(I,k); it; ++it){
        double v = it.value();
        int i = it.row();   // row index
        int j = it.col();   // col index (here it is equal to k)
        long long int p,q;
        if(dbl2frac(1.0e-12,v,&p,&q)){
          printf("## THE TIME AS COME TO REPLACE DBL2FRAC CF tools/gen_lag2bezfuncs.cxx\n");
          printf("Lowering METRIS_MAX_DEG_JACOBIAN will work if %d is not necessary\n",METRIS_MAX_DEG_JACOBIAN);
          printf("Failure at deg  %d on coeff = %30.23f \n",ideg,v);
          printf("Call dbgl2frac verbose:\n");
          dbl2frac(1.0e-12,v,&p,&q,true);
          exit(1);
        }
        if(abs((double ) p/ (double) q - v) > 1.0e-6){
          printf("## LARGE ERROR %15.7e %15.7e\n",(double ) p/ (double) q,v);
        }
        pcoefs[i].push_back(p);
        qcoefs[i].push_back(q);
        idxs[i].push_back(j);
      }
    }
    #else
    for (int ii = 0; ii < npp; ++ii){
      for (int jj = 0; jj < npp; jj++){
        double v = invM(ii,jj);
        long long int p,q;
        if(dbl2frac(1.0e-12,v,&p,&q)){
          printf("## THE TIME AS COME TO REPLACE DBL2FRAC CF tools/gen_lag2bezfuncs.cxx\n");
          printf("Lowering METRIS_MAX_DEG_JACOBIAN will work if %d is not necessary\n",METRIS_MAX_DEG_JACOBIAN);
          printf("Failure at deg  %d on coeff = %30.23f \n",ideg,v);
          printf("Call dbgl2frac verbose:\n");
          dbl2frac(1.0e-12,v,&p,&q,true);
          exit(1);
        }
        if(abs((double ) p/ (double) q - v) > 1.0e-6){
          printf("## LARGE ERROR %15.7e %15.7e\n",(double ) p/ (double) q,v);
        }
        pcoefs[ii].push_back(p);
        qcoefs[ii].push_back(q);
        idxs[ii].push_back(jj);
      }
    }
    #endif

    for(int isz = 0; isz < nszfld; isz++){
      int szfld = szflds[isz];

      str << "template<> void lag2bez"<<tdim<<"<"<<ideg<<","<<szfld<<">"<<
                               "(const int* __restrict__ lfld,\n\
                                const dblAr2& __restrict__ rfld0,\n\
                                dblAr2& __restrict__ rfld1){\n";

      if(szfld > 1){
        str << "  for(int i=0; i< "<<szfld<<"; i++){\n";
      }
      for(int i=0; i < tdim + 1 ;i++){
        snprintf(i_s,4,"%3d",i);
        if(szfld > 1){
          str << "    rfld1[lfld["<<i_s<<"]][i] = rfld0[lfld["<<i_s<<"]][i];\n";
        }else{
          str << "  rfld1[lfld["<<i_s<<"]][0] = rfld0[lfld["<<i_s<<"]][0];\n";
        }
      }
  
      for(int i = tdim + 1; i < npp; i++){
        snprintf(i_s,4,"%3d",i);
          //      double v = coeffs[i][0];
        int j = idxs[i][0];
        long long int p = pcoefs[i][0];
        long long int q = qcoefs[i][0];
          //      snprintf(v_s,64,"%22.16f",v);
          //      snprintf(p_s,64,"%22.16f",p);
          //      snprintf(q_s,64,"%22.16f",q);
        snprintf(j_s,4,"%3d",j);
  
        if(szfld > 1){
          //        str << "    rfld1[lfld["<<i_s<<"]][i] = "<<v_s<<"*rfld0[lfld["<<j_s<<"]][i]";
          str << "    rfld1[lfld["<<i_s<<"]][i] = "<<p<<"*rfld0[lfld["<<j_s<<"]][i]/"<<q;
          for(int j = 1; j < pcoefs[i].size(); j++){
          //          v = coeffs[i][j];
          //          snprintf(v_s,64,"%22.16f",v);
          //          snprintf(p_s,64,"%22.16f",v);
          //          snprintf(q_s,64,"%22.16f",v);
            p = pcoefs[i][j];
            q = qcoefs[i][j];
            int jj = idxs[i][j];
            snprintf(j_s,4,"%3d",jj);
          //          str<<" + "<<v_s<<"*rfld0[lfld["<<j_s<<"]][i]";
            str<<" + "<<p<<"*rfld0[lfld["<<j_s<<"]][i]/"<<q;
          }
        }else{
          //        str << "  rfld1[lfld["<<i_s<<"]][0] = "<<v_s<<"*rfld0[lfld["<<j_s<<"]][0]";
          str << "  rfld1[lfld["<<i_s<<"]][0] = "<<p<<"*rfld0[lfld["<<j_s<<"]][0]/"<<q;
          //        if(iprt> 0)printf("Debug within loop i = %d size %d \n",i,(int)coeffs[i].size());
          for(int j = 1; j < pcoefs[i].size(); j++){
          //          if(iprt > 0)printf("%d / %d \n",j,(int)(coeffs[i].size()-1));
          //          v = coeffs[i][j];
          //          snprintf(v_s,64,"%22.16f",v);
            p = pcoefs[i][j];
            q = qcoefs[i][j];
            int jj = idxs[i][j];
            snprintf(j_s,4,"%3d",jj);
          //          if(iprt> 0)printf("idx %d \n",idxs[i][j]);
          //          str<<" + "<<v_s<<"*rfld0[lfld["<<j_s<<"]][0]";
            str<<" + "<<p<<"*rfld0[lfld["<<j_s<<"]][0]/"<<q;
          }
        }
      	str<<";\n";
      }
      if(szfld > 1){
        str << "  }\n";
      }
      str<<"}\n\n";
    }
  }CT_FOR1(ideg);
  //  str << "  }\n"; // close constructor
  //  str << "  ALWAYS_INLINE double* operator[](int i)const {\n";
  //  str << "    return arr[i];\n";
  //  str << "  }\n";
  ////  str << "  const double arr[tetnpps[ideg]][tetnpps[ideg]];";
  //  str << "  dblAr2 arr;\n";
  //  str<<"};\n"; // close struct
  //  str << "#endif";


  //  std::ofstream f;
  //  char fname[64]; snprintf(fname,64,"codegen_lag2bez%d.cxx",tdim);
  //  f.open(fname, std::ios::out);
  //  f << str.str();
  //  f.close();
}

template void gen_lag2bez<1>(std::ostringstream &str);
template void gen_lag2bez<2>(std::ostringstream &str);
template void gen_lag2bez<3>(std::ostringstream &str);

int main(int argc, char **argv){
	char fname[64];
  std::ofstream f;

  std::ostringstream strh;
  strh<<"#ifndef __LAG2BEZ__\n";
  strh<<"#define __LAG2BEZ__\n\n";

//  strh<<"#include \"common_includes.hxx\"\n\n";
//  strh<<"#include \"ho_constants.hxx\"\n";

  strh << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  strh << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  strh << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  strh << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";
  strh<<"#include \"types.hxx\"\n\n";
  strh<<"namespace Metris{\n\n";

  for(int tdim = 1; tdim <= 3; tdim++){
  	strh << "template<int ideg, int szfld>\n";
  	strh << "void lag2bez"<<tdim<<
  	          "(const int* __restrict__ lfld,\n\
 								const dblAr2& __restrict__ rfld0,\n\
 								dblAr2& __restrict__ rfld1);\n";
  }
  strh<<"\n} // End namespace\n\n";
  strh<<"\n#endif\n";
  snprintf(fname,64,"codegen_lag2bez.hxx");
  f.open(fname, std::ios::out);
  f << strh.str();
  f.close();




  std::ostringstream strc;

  strc << "//Metris: high-order metric-based non-manifold tetrahedral remesher\n";
  strc << "//Copyright (C) 2023-2024, Massachusetts Institute of Technology\n";
  strc << "//Licensed under The GNU Lesser General Public License, version 2.1\n";
  strc << "//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php\n\n";
  strc<<"#include \"codegen_lag2bez.hxx\"\n\n";
  strc<<"namespace Metris{\n\n";
  gen_lag2bez<1>(strc);
  gen_lag2bez<2>(strc);
  gen_lag2bez<3>(strc);
  strc<<"\n} // End namespace\n\n";
  snprintf(fname,64,"codegen_lag2bez.cxx");
  f.open(fname, std::ios::out);
  f << strc.str();
  f.close();
  return 0;
}












