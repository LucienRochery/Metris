//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __MSH_DEGELEV__
#define __MSH_DEGELEV__


#include "ho_constants.hxx"
#include "low_eval.hxx"
#include "Mesh/MeshFwd.hxx"
#include "Mesh/MeshBase.hxx"


namespace Metris{

// tdim: topo dimension
// gdim: geometric dimension
template<int mshdeg, int tardeg, int tdim, int gdim>
struct gen_nodes_degelev{
  gen_nodes_degelev() = delete;
  gen_nodes_degelev(const MeshBase &msh, int ientt, double* newpt){

    static_assert((tdim <= 3 && "## gen_nodes_degelev TOPO DIM > 3 NOT IMPLEMENTED \n"));

    constexpr auto entnpps = ENTNPPS(tdim);
    constexpr auto ordent  = ORDELT(tdim);
    constexpr int npptar = entnpps[tardeg];
    constexpr int nppsrc = entnpps[mshdeg];

    constexpr auto eval = tdim == 1 ? eval1<gdim,mshdeg> 
                        : tdim == 2 ? eval2<gdim,mshdeg> 
                                    : eval3<gdim,mshdeg> ;


    double bary[tdim+1];

    if (msh.getBasis() == FEBasis::Lagrange){
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      // const int ntmp = tetnpps[tardeg];
      for(int irnk = 0; irnk < npptar; irnk++){
        for(int i = 0; i < tdim+1 ;i++) {
          bary[i] = ordent[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval(msh.coord,&ent2poi(ientt,0),FEBasis::Lagrange,
             DifVar::None,DifVar::None,bary,&newpt[gdim*irnk],NULL,NULL);
      }
      return;
    }else{
      for(int irnk = 0; irnk < npptar; irnk++){
        
        for(int i = 0;i<gdim;i++){
          newpt[gdim*irnk+i] = 0.0;
        }
      }
      // -- Loop over original degree control points
      for(int irnk1 = 0; irnk1 < nppsrc ;irnk1++){
        int ii[tdim+1],jj[tdim+1];
        for(int i = 0; i < tdim + 1; i++){
          ii[i] = ordent[mshdeg][irnk1][i];
        }
        //if constexpr(tdim == 1){
        //  ii[0] = ordedg.s[mshdeg][irnk1][0];
        //  ii[1] = ordedg.s[mshdeg][irnk1][1];
        //}else if(tdim == 2){
        //  ii[0] = ordfac.s[mshdeg][irnk1][0];
        //  ii[1] = ordfac.s[mshdeg][irnk1][1];
        //  ii[2] = ordfac.s[mshdeg][irnk1][2];
        //}else if(tdim == 3){
        //  ii[0] = ordtet.s[mshdeg][irnk1][0];
        //  ii[1] = ordtet.s[mshdeg][irnk1][1];
        //  ii[2] = ordtet.s[mshdeg][irnk1][2];
        //  ii[3] = ordtet.s[mshdeg][irnk1][3];
        //}
        // -- Multiply by the sum of all degree tardeg - msh.curdeg Bernsteins
        // -- Distribute to a + b with coeff, this is in factor of P_a
        for(int irnk2 = 0; irnk2<entnpps[tardeg-msh.curdeg];irnk2++){
          int irnk3;
          double coef;

          for(int i = 0; i < tdim + 1; i++){
            jj[i] = ordent[tardeg-mshdeg][irnk2][i];
          }

          if constexpr(tdim == 1){
            //jj[0] = ordedg.s[tardeg-mshdeg][irnk2][0];
            //jj[1] = ordedg.s[tardeg-mshdeg][irnk2][1];
            irnk3 = mul2nod(ii[0]+jj[0],ii[1]+jj[1]);
            coef = cbzedg.s[tardeg][irnk3] / 
                     (cbzedg.s[mshdeg][irnk1]*cbzedg.s[tardeg-mshdeg][irnk2]);
            
            for(int i = 0;i<gdim;i++){
              newpt[gdim*irnk3+i] 
               += msh.coord[msh.edg2poi(ientt,irnk1)][i] / coef;
            }
          }else if(tdim == 2){
            //jj[0] = ordfac.s[tardeg-mshdeg][irnk2][0];
            //jj[1] = ordfac.s[tardeg-mshdeg][irnk2][1];
            //jj[2] = ordfac.s[tardeg-mshdeg][irnk2][2];
            irnk3 = mul2nod(ii[0]+jj[0],ii[1]+jj[1],ii[2]+jj[2]);
            coef = cbzfac.s[tardeg][irnk3] / 
                     (cbzfac.s[mshdeg][irnk1]*cbzfac.s[tardeg-mshdeg][irnk2]);
            
            for(int i = 0;i<gdim;i++){
              newpt[gdim*irnk3+i] 
               += msh.coord[msh.fac2poi(ientt,irnk1)][i] / coef;
            }
          }else if(tdim == 3){
            //jj[0] = ordtet.s[tardeg-mshdeg][irnk2][0];
            //jj[1] = ordtet.s[tardeg-mshdeg][irnk2][1];
            //jj[2] = ordtet.s[tardeg-mshdeg][irnk2][2];
            //jj[3] = ordtet.s[tardeg-mshdeg][irnk2][3];
            irnk3 = mul2nod(ii[0]+jj[0],ii[1]+jj[1],ii[2]+jj[2],ii[3]+jj[3]);
            coef = cbztet.s[tardeg][irnk3] / 
                     (cbztet.s[mshdeg][irnk1]*cbztet.s[tardeg-mshdeg][irnk2]);
            
            for(int i = 0;i<gdim;i++){
              newpt[gdim*irnk3+i] 
               += msh.coord[msh.tet2poi(ientt,irnk1)][i] / coef;
            }
          }
        }
      }
    } 

  }
};

// tdim: topo dimension
// gdim: geometric dimension; really just the stride of newpt
// lbpoi: the ibpois of the points; this may require some linked list walking 
// This is an input as it will be reused by the called. 
template<int mshdeg, int tardeg, int tdim, int uvdim>
struct gen_uvs_degelev{
  gen_uvs_degelev() = delete;
  gen_uvs_degelev(const MeshBase &msh, int* lbpoi, double* newuv){

    static_assert((tdim <= 2 && "gen_uvs_degelev TOPO DIM > 2 NOT IMPLEMENTED"));

    constexpr int npptar = (tdim == 1) ? edgnpps[tardeg] : facnpps[tardeg];

    double bary[tdim+1];
    for(int irnk = 0; irnk < npptar; irnk++){
      if constexpr(tdim == 1){
        for(int i = 0; i < tdim + 1 ;i++) {
          bary[i] = ordedg.s[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval1<1,mshdeg>(msh.bpo2rbi,lbpoi,FEBasis::Lagrange,DifVar::None,DifVar::None,bary,&newuv[uvdim*irnk],NULL,NULL);
      }else {
        for(int i = 0; i < tdim + 1 ;i++) {
          bary[i] = ordfac.s[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval2<2,mshdeg>(msh.bpo2rbi,lbpoi,FEBasis::Lagrange,DifVar::None,DifVar::None,bary,&newuv[uvdim*irnk],NULL,NULL);
      }
    }

  }
};


/* Deprecated, now use met.genMetNodes();
template<int mshdeg, int tardeg, int tdim, int gdim>
struct gen_met_degelev{
  gen_met_degelev() = delete;
  gen_met_degelev(const MeshBase &msh, MetricFieldFE &met, int ientt, double* newmet){

    static_assert((tdim <= 3 && "## gen_nodes_degelev TOPO DIM > 3 NOT IMPLEMENTED \n"));

    METRIS_ENFORCE_MSG(met.getBasis() == FEBasis::Lagrange, "Add lag2bez if intending to use this with Bézier metric");

    constexpr int npptar = (tdim == 1) ? edgnpps[tardeg]
                         :((tdim == 2) ? facnpps[tardeg]
                                       : tetnpps[tardeg]);
    constexpr int nppsrc = (tdim == 1) ? edgnpps[mshdeg]
                         :((tdim == 2) ? facnpps[mshdeg]
                                       : tetnpps[mshdeg]);

    constexpr int nnmet = (gdim*(gdim+1))/2;

    double bary[tdim+1];

    for(int irnk = 0; irnk < npptar; irnk++){
      if constexpr(tdim == 1){
        for(int i = 0; i < 2 ;i++) {
          bary[i] = ordedg.s[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval1<nnmet,mshdeg>(msh.met,&msh.edg2poi(ientt,0),met.getBasis(),0,0,bary,&newmet[nnmet*irnk],NULL,NULL);
      }else if (tdim == 2){
        for(int i = 0; i < 3 ;i++) {
          bary[i] = ordfac.s[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval2<nnmet,mshdeg>(msh.met,&msh.fac2poi(ientt,0),met.getBasis(),0,0,bary,&newmet[nnmet*irnk],NULL,NULL);
      }else if (tdim == 3){
        for(int i = 0; i < 4 ;i++) {
          bary[i] = ordtet.s[tardeg][irnk][i] / (1.0*tardeg);
        }
        eval3<nnmet,mshdeg>(msh.met,&msh.tet2poi(ientt,0),met.getBasis(),0,0,bary,&newmet[nnmet*irnk],NULL,NULL);
      }
    }

  }
};
*/



/*
Routine reworked by mid-June after surface topological handling was rewritten. 
It now seems more straightforward to proceed from edges, to triangles, to tetrahedra rather than the opposite. 
Let us also save on evaluations by first checking the neighbours and only then generating as needed. 
This implements topological and (u,v) guess init of boundary points.
*/
template<class MetricFieldType, int mshdeg, int tardeg>
void deg_elevate(Mesh<MetricFieldType> &msh);


} // End namespace

#endif
