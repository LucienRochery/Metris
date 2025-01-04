//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/utils.hxx"


namespace Metris{

template<int gdim, int nvar>
void getmet_dbl2SurS(const double*__restrict__ met, const double*__restrict__ dmet, 
									SANS::SurrealS<nvar,double>*__restrict__ metS){
	constexpr int nnmet = (gdim*(gdim+1))/2;
	for(int ii = 0; ii < nnmet ;ii++){
		metS[ii].value() = met[ii];
		if(dmet == NULL) continue;
		for(int jj = 0; jj < nvar; jj++){
			metS[ii].deriv(jj) = dmet[nnmet*jj + ii];
		}
	}
}

template void getmet_dbl2SurS<2,2>(const double*__restrict__ met, const double*__restrict__ dmet, 
									SANS::SurrealS<2,double>*__restrict__ metS);
template void getmet_dbl2SurS<3,2>(const double*__restrict__ met, const double*__restrict__ dmet, 
									SANS::SurrealS<2,double>*__restrict__ metS);
template void getmet_dbl2SurS<3,3>(const double*__restrict__ met, const double*__restrict__ dmet, 
                  SANS::SurrealS<3,double>*__restrict__ metS);



template<int gdim, int nvar>
void getmet_SurS2dbl(const SANS::SurrealS<nvar,double>*__restrict__ metS,
										 double*__restrict__ met, double*__restrict__ dmet ){
	constexpr int nnmet = (gdim*(gdim+1))/2;
	for(int ii = 0; ii < nnmet ;ii++){
		met[ii] = metS[ii].value();
		if(dmet == NULL) continue;
		for(int jj = 0; jj < nvar; jj++){
			dmet[nnmet*jj + ii] = metS[ii].deriv(jj);
		}
	}
}

template void getmet_SurS2dbl<2,2>(const SANS::SurrealS<2,double>*__restrict__ metS,
										 double*__restrict__ met, double*__restrict__ dmet);
template void getmet_SurS2dbl<3,2>(const SANS::SurrealS<2,double>*__restrict__ metS,
										 double*__restrict__ met, double*__restrict__ dmet);
template void getmet_SurS2dbl<3,3>(const SANS::SurrealS<3,double>*__restrict__ metS,
                     double*__restrict__ met, double*__restrict__ dmet);
}// End namsepace
