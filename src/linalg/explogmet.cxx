//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/explogmet.hxx"
#include "../linalg/eigen.hxx"
#include "../aux_exceptions.hxx"

#include "../low_geo.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/utils.hxx"


namespace Metris{


// -----------------------------------------------------------------------------
// Write log(met) in lmet
template<int ndim, typename T>
void getlogmet_cpy(const T* __restrict__ met, T* __restrict__ lmet){
	static_assert(ndim == 2 || ndim == 3);

	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);
	if(eigval[0] < 1.0e-16) METRIS_THROW_MSG(RealExcept(),"Negative eigenvalues min = "<<
			eigval[0]<<"\n");

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = log(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,lmet);
}

// -----------------------------------------------------------------------------
// Write log(met) as sum to lmet already initialized
template<int ndim, typename T>
void getlogmet_sum(const T* __restrict__ met, T* __restrict__ lmet){
	static_assert(ndim == 2 || ndim == 3);

	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);
	if(eigval[0] < 1.0e-16) METRIS_THROW_MSG(RealExcept(),"Negative eigenvalues min = "<<
			eigval[0]<<"\n");

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = log(eigval[ii]);

	eig2met_sum<ndim,T>(eigval,eigvec,lmet);
}


// -----------------------------------------------------------------------------
// Replace met with lmet
template<int ndim, typename T>
void getlogmet_inp(T *met){
	static_assert(ndim == 2 || ndim == 3);

	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);
	if(eigval[0] < 1.0e-16){
    printf("Invalid metric: \n");
    int nnmet = (ndim*(ndim+1))/2;
    if constexpr (std::is_same<T, double>::value){
      for(int ii = 0 ; ii < nnmet; ii++) printf(" %23.15e ",met[ii]); 
    }else{
      for(int ii = 0 ; ii < nnmet; ii++) std::cout<<met[ii]<<" ";
    }
    std::cout<<"\n";

    if constexpr(std::is_same<T,double>::value){
      std::cout<<"eigvals:";
      dblAr1(ndim,eigval).print();
    }

    METRIS_THROW_MSG(RealExcept(),"Negative eigenvalues");
  }

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = log(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,met);
}


template <int ndim, typename T>
void getexpmet_inp(T* met){
	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = exp(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,met);
}


template void getlogmet_cpy<2,double>(const double* __restrict__, double *__restrict__);
template void getlogmet_cpy<3,double>(const double* __restrict__, double *__restrict__);

template void getlogmet_sum<2,double>(const double* __restrict__, double *__restrict__);
template void getlogmet_sum<3,double>(const double* __restrict__, double *__restrict__);

template void getlogmet_inp<2,double>(double*);
template void getlogmet_inp<2,SANS::SurrealS<2,double>>(SANS::SurrealS<2,double>*);
template void getlogmet_inp<3,double>(double*);
template void getlogmet_inp<3,SANS::SurrealS<3,double>>(SANS::SurrealS<3,double>*);

template void getexpmet_inp<2,double>(double*);
template void getexpmet_inp<2,SANS::SurrealS<2,double>>(SANS::SurrealS<2,double>*);
template void getexpmet_inp<3,double>(double*);
template void getexpmet_inp<3,SANS::SurrealS<3,double>>(SANS::SurrealS<3,double>*);




template <int n>
void getexpmet_cpy(const double* met ,double*  __restrict__ expm, 
                   double tol, int iscal){
  constexpr int nnmet = (n*(n+1))/2;
  double iterm[2][nnmet];
  int iwhich,inmet;

  double nrm = getnrml2<nnmet>(met);
  
  int powscal = 0;
  double scalfac = 1;
  while(nrm > 0.25 && iscal){
    powscal++;
    scalfac /= 2;
    nrm /= 4;
  }
  for(inmet = 0; inmet < nnmet; inmet++){
    iterm[0][inmet] = met[inmet]*scalfac;
  }
	//     --- I + iterm
  for(int ii = 0; ii < n; ii++){
    for(int jj = ii ; jj < n; jj++){
      expm[sym2idx(ii,jj)] = (ii == jj) + iterm[0][sym2idx(ii,jj)];
    }
  }
  //expm[0] = 1 + iterm[0][0];
  //expm[1] =     iterm[0][1];
  //expm[2] = 1 + iterm[0][2];
  //expm[3] =     iterm[0][3];
  //expm[4] =     iterm[0][4];
  //expm[5] = 1 + iterm[0][5];
 
  iwhich = 0;
	int niter = 1;
  while(nrm > tol*tol && niter < 1000){
    niter++;

    symXsymsub_fac<n>(iterm[iwhich],met,scalfac/niter,iterm[1-iwhich]);

    iwhich = 1-iwhich;
    
    for(inmet = 0; inmet < nnmet; inmet++){
      expm[inmet] += iterm[iwhich][inmet];
    }
    
    nrm = geterrl2<nnmet>(iterm[0],iterm[1]);
  }

  if(nrm > tol*tol) METRIS_THROW(RealExcept());

  if(powscal == 0) return;
  for(inmet = 0; inmet < nnmet; inmet++){
    iterm[iwhich][inmet] = expm[inmet];
  }
  for(niter = 0; niter < powscal; niter ++){
	//   --- T^(2n)
    symXsymsub<n>(iterm[iwhich],iterm[iwhich],iterm[1-iwhich]);

    iwhich = 1-iwhich;
  }

  for(inmet = 0; inmet < nnmet; inmet++){
    expm[inmet] = iterm[iwhich][inmet];
  }

}

template void getexpmet_cpy<2>(const double* met ,double*  __restrict__ expm, 
                   double tol, int iscal);
template void getexpmet_cpy<3>(const double* met ,double*  __restrict__ expm, 
                   double tol, int iscal);


void getexpmet_cpy_d(const double* met ,const double* dmet, 
                    double*  __restrict__ expm,double*  __restrict__ dexp, 
                    double tol, int iscal){
	double nrm,scalfac, iterm[2][6],gterm[2][3][6];
	int niter,powscal,iwhich,igrad,inmet;
	niter = 1;
  nrm = getnrml2<6>(met);
  
  powscal = 0;
  scalfac = 1;
  while(nrm > 0.25 && iscal){
   powscal++;
   scalfac /= 2;
   nrm /= 4;
  }
  for(inmet = 0; inmet < 6; inmet++){
    iterm[0] [inmet] = met[inmet]*scalfac;
    dexp [0*6+inmet] = gterm[0][0][inmet] = dmet[0*6+inmet]*scalfac;
    dexp [1*6+inmet] = gterm[0][1][inmet] = dmet[1*6+inmet]*scalfac;
    dexp [2*6+inmet] = gterm[0][2][inmet] = dmet[2*6+inmet]*scalfac;
  }
	//     --- I + iterm
  expm[0] = 1 + iterm[0][0];
  expm[1] =     iterm[0][1];
  expm[2] = 1 + iterm[0][2];
  expm[3] =     iterm[0][3];
  expm[4] =     iterm[0][4];
  expm[5] = 1 + iterm[0][5];

  iwhich = 0;
  while(nrm > tol*tol && niter < 1000){
   niter++;

    symXsymsub_fac<3>(iterm[iwhich],met,scalfac/niter,iterm[1-iwhich]);
    for(igrad = 0; igrad < 3; igrad++){
      symXsymsub_fac<3>(iterm[iwhich],&dmet[6*igrad],scalfac/niter,gterm[1-iwhich][igrad]);
      symXsymadd_fac<3>(gterm[iwhich][igrad],met,scalfac/niter,gterm[1-iwhich][igrad]);
    }

    iwhich = 1-iwhich;

   for(inmet = 0; inmet < 6; inmet++){
     expm    [inmet] += iterm[iwhich]   [inmet];
      dexp[0*6+inmet] += gterm[iwhich][0][inmet];
      dexp[1*6+inmet] += gterm[iwhich][1][inmet];
      dexp[2*6+inmet] += gterm[iwhich][2][inmet];
   }

    nrm = geterrl2<6>(iterm[0],iterm[1]);

    //for(i = 0; i < 3; i++){
    //  double tmp = geterrl2<6>(gterm[0][i],gterm[1][i]);
    //  nrm = nrm > tmp ? nrm : tmp;
    //}
  }
  if(nrm > tol*tol) METRIS_THROW(RealExcept());

  if(powscal == 0) return;
  for(inmet = 0; inmet < 6; inmet++){
    iterm[iwhich]   [inmet] = expm   [inmet];
    gterm[iwhich][0][inmet] = dexp[0*6+inmet];
    gterm[iwhich][1][inmet] = dexp[1*6+inmet];
    gterm[iwhich][2][inmet] = dexp[2*6+inmet];
  }

  for(niter = 0; niter < powscal; niter ++){
		//   --- T^(2n)
    symXsymsub<3>(iterm[iwhich],iterm[iwhich],iterm[1-iwhich]);

		//   --- DT^(2n)
    ABpBA3symsub(iterm[iwhich],gterm[iwhich][0],gterm[1-iwhich][0]);
    ABpBA3symsub(iterm[iwhich],gterm[iwhich][1],gterm[1-iwhich][1]);
    ABpBA3symsub(iterm[iwhich],gterm[iwhich][2],gterm[1-iwhich][2]);

    iwhich = 1-iwhich;
  }

  for(inmet = 0; inmet < 6; inmet++){
    expm    [inmet] = iterm[iwhich]   [inmet];
    dexp[0*6+inmet] = gterm[iwhich][0][inmet];
    dexp[1*6+inmet] = gterm[iwhich][1][inmet];
    dexp[2*6+inmet] = gterm[iwhich][2][inmet];
  }

}



// Using surreals. For validation purposes. 10% slower than direct approach.
void getexpmet_cpy_dS(const double* met ,const double* dmet,  
                     double*  __restrict__ expm,double*  __restrict__ dexp, 
                     double tol, int iscal){
	double nrm,scalfac;
	SANS::SurrealS<3,double> iterm[2][6];
	SANS::SurrealS<3,double> metS[6];
  getmet_dbl2SurS<3>(met,dmet,metS);

	//for(int i = 0 ; i < 6 ; i++){
	//	metS[i].value() = met[i];
	//	for(int j = 0; j < 3 ; j++){
	//		metS[i].deriv(j) = dmet[j*6 + i];
	//	}
	//}

	SANS::SurrealS<3,double> expmS[6];
	int niter,powscal,iwhich,inmet;
	niter = 1;
  nrm = getnrml2<6>(met);
  
  powscal = 0;
  scalfac = 1;
  while(nrm > 0.25 && iscal){
   powscal++;
   scalfac /= 2;
   nrm /= 4;
  }
  for(inmet = 0; inmet < 6; inmet++){
    iterm[0][inmet] = metS[inmet]*scalfac;
  }

	//     --- I + iterm
  expmS[0] = iterm[0][0] + 1 ;
  expmS[1] = iterm[0][1]     ;
  expmS[2] = iterm[0][2] + 1 ;
  expmS[3] = iterm[0][3]     ;
  expmS[4] = iterm[0][4]     ;
  expmS[5] = iterm[0][5] + 1 ;


  iwhich = 0;
  while(nrm > tol*tol && niter < 1000){
   niter++;

    symXsymsub_fac<3>(iterm[iwhich],metS,scalfac/niter,iterm[1-iwhich]);

    iwhich = 1-iwhich;

   for(inmet = 0; inmet < 6; inmet++){
     expmS[inmet] += iterm[iwhich][inmet];
   }

    nrm = geterrl2<6,3>(iterm[0],iterm[1]);

  }
  if(nrm > tol*tol) METRIS_THROW(RealExcept());

  if(powscal == 0){
   for(inmet = 0; inmet < 6; inmet++){
     expm    [inmet] = expmS[inmet].value();
     dexp[0*6+inmet] = expmS[inmet].deriv(0);
     dexp[1*6+inmet] = expmS[inmet].deriv(1);
     dexp[2*6+inmet] = expmS[inmet].deriv(2);
   }
   return;
  }
  for(inmet = 0; inmet < 6; inmet++){
    iterm[iwhich][inmet] = expmS[inmet];
  }

  for(niter = 0; niter < powscal; niter ++){
		//   --- T^(2n)
    symXsymsub<3>(iterm[iwhich],iterm[iwhich],iterm[1-iwhich]);
    iwhich = 1-iwhich;
  }



  for(inmet = 0; inmet < 6; inmet++){
    expm    [inmet] = iterm[iwhich][inmet].value();
    dexp[0*6+inmet] = iterm[iwhich][inmet].deriv(0);
    dexp[1*6+inmet] = iterm[iwhich][inmet].deriv(1);
    dexp[2*6+inmet] = iterm[iwhich][inmet].deriv(2);
  }


}


} // End namespace

