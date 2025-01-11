//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Mesh/Mesh.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../metris_constants.hxx"
#include "../ho_constants.hxx"
#include "../aux_utils.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_topo.hxx"
#include "../low_lenedg.hxx"
#include "../low_geo.hxx"
#include "../low_ccoef.hxx"

#include "../Localization/msh_localization.hxx"


namespace Metris{

template<class MFT>
int Mesh<MFT>::newpoitopo(int tdimn, int ientt){
  int ipoin = MeshBase::newpoitopo(tdimn, ientt);

  // No need to add guesses here, that'll be done when calling interpMetBack
  for(int tdim2 = 1; tdim2 <= this->get_tdim(); tdim2++) 
    poi2bak(ipoin,tdim2-1) = -1;

  return ipoin;
}


template<>
void Mesh<MetricFieldAnalytical>::initialize(MetrisAPI *data, MeshBack &bak, 
  #ifdef NDEBUG 
  const 
  #endif
  MetrisParameters &param){
  this->initializeCommon(data,bak,param);

  if(param.ianamet >= 0){
    met.setAnalyticalMetric(param.ianamet);
  }else{
    METRIS_ASSERT(param.anamet_ptr != NULL);
    met.setAnalyticalMetric(param.anamet_ptr);
  }
  if(param.scaleMet){
    if(param.iverb >= 1) 
      printf("-- Front scaling metric by %15.7e\n", param.metScale);
    met.normalize(param.metScale);
  }

  met.forceBasisFlag(FEBasis::Lagrange);
  met.forceSpaceFlag(MetSpace::Exp);

  FEBasis mshbas0 = this->ibasis;
  this->setBasis(FEBasis::Lagrange);
  for(int ipoin = 0; ipoin < npoin; ipoin++){
    met.getMetPhys(AsDeg::Pk,DifVar::None,MetSpace::Exp,NULL,
                   coord[ipoin],met[ipoin],NULL);
  }
  this->setBasis(mshbas0);

}



template<>
void Mesh<MetricFieldFE>::initialize(MetrisAPI *data, MeshBack &bak, 
  #ifdef NDEBUG 
  const 
  #endif
  MetrisParameters &param){
  this->initializeCommon(data,bak,param);

  poi2bak.fill(-1);

  if(data == NULL && !param.inpBack){ // Copy case

    met = bak.met;

    for(int tdim = 1; tdim <= get_tdim(); tdim++){
      int nentt = this->nentt(tdim); 
      const intAr2 &ent2poi = this->ent2poi(tdim);
      int nnode = this->nnode(tdim);

      for(int ii = 0 ; ii < nentt; ii++){
        if(isdeadent(ii,ent2poi)) continue;
        for(int jj = 0; jj < nnode; jj++){
          int ip = ent2poi(ii,jj);
          poi2bak(ip,tdim-1) = ii;
        }
      }
    }


  }else{

    if(param.iverb >= 1) std::cout<<"-- Back metric interpolation\n";
    CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bak.curdeg == bdeg){
      interpFrontBack<MetricFieldFE,bdeg>(*this,bak);
    }}CT_FOR1(bdeg);

  }


}


template<class MFT>
void Mesh<MFT>::initializeCommon(MetrisAPI *data, MeshBack &bak, 
  #ifdef NDEBUG 
  const 
  #endif
  MetrisParameters &param){
  this->param = &param;
  this->bak = &bak;

  // The back mesh has curdeg = strdeg regardless of target degree 
  this->strdeg = MAX(param.usrTarDeg,bak.curdeg);

  if(data == NULL && !param.inpBack){
    // In this case, simply copy from back mesh 
    MeshBase &dum = *this;
    dum = (MeshBase &) bak; // MeshBase::operator= 
  }else{
    MeshBase::initialize(data,param);
  }
  
  //this->cfa2tag.allocate(METRIS_MAXTAGS, this->ncadfa, true);
  //this->ced2tag.allocate(METRIS_MAXTAGS, this->ncaded, true);
  //this->cno2tag.allocate(METRIS_MAXTAGS, this->ncadno, true);

  //this->cfa2tag.fill(METRIS_MAXTAGS, this->ncadfa,0);
  //this->ced2tag.fill(METRIS_MAXTAGS, this->ncaded,0);
  //this->cno2tag.fill(METRIS_MAXTAGS, this->ncadno,0);

}

template<class MFT>
void Mesh_getEnttMemCosts(const Mesh<MFT> &msh, 
                        int *memCostPpoi, int *memCostPbpo, int *memCostPedg, 
                        int *memCostPfac, int *memCostPelt){
  msh.MeshBase::getEnttMemCosts(memCostPpoi,memCostPbpo,memCostPedg,memCostPfac,memCostPelt);

  int memCostPdbl = sizeof(double);
  int memCostPint = sizeof(int);
  
  int nnmet = (msh.idim*(msh.idim + 1))/2;

  *memCostPpoi += nnmet*memCostPdbl  /* met     */
                +     1*memCostPint; /* poi2bak */
}


template
void Mesh_getEnttMemCosts<MetricFieldFE>(const Mesh<MetricFieldFE>&,
                                      int*,int*,int*,int*,int*);
template
void Mesh_getEnttMemCosts<MetricFieldAnalytical>(const Mesh<MetricFieldAnalytical>&,
                                      int*,int*,int*,int*,int*);




template<class MFT>
void Mesh<MFT>::set_npoin(int npoin, bool skipallocf){
  MeshMetric<MFT>::set_npoin(npoin, skipallocf);
  // Only allocate needed dimensions. Say dim = 2, allocate for edges and tri. 
  poi2bak.allocate(this->mpoin, this->idim);
  poi2bak.set_n(this->npoin); 
}


template<class MFT>
void Mesh<MFT>::set_nentt(int tdim, int nentt, bool skipallocf){
  switch(tdim){
  case(-1):
    MeshBase::set_nbpoi(nentt);
    break;
  case(0):
    set_npoin(nentt,skipallocf);
    break;
  case(1):
    MeshBase::set_nedge(nentt,skipallocf);
    break;
  case(2):
    MeshBase::set_nface(nentt,skipallocf);
    break;
  case(3):
    MeshBase::set_nelem(nentt,skipallocf);
    break;
  }
}


// tdim is the dimension of iseed (front elt), not of the point !
// though in first implementation, these should be the same, i.e. lowest. 
// iref and algnd only necessary if doing boundary localization. 
template<class MFT>
int Mesh<MFT>::interpMetBack(int ipoin, int tdim, int iseed, 
                             int iref, const double* algnd){
  METRIS_ASSERT_MSG(tdim == this->getpoitdim(ipoin) || this->getpoitdim(ipoin) == 0,
    "seed is dim "<<tdim<<" point is "<<this->getpoitdim(ipoin)
    << " ipoin = "<<ipoin );
  int pdim = this->getpoitdim(ipoin);

  if(this->param->iverb >= 3) 
    printf("\n  -- START interpMetBack ipoin = %d iseed = %d tdim %d \n",ipoin,
      iseed,tdim);

  METRIS_ASSERT(tdim == this->idim || (algnd != NULL && iref >= 0));

  METRIS_ASSERT_MSG(ipoin >= 0 && ipoin < this->npoin, 
    "interpMetBack ipoin out of bounds "<<ipoin<<" < ? "<<this->npoin);

  int ierro = this->interpMetBack0(ipoin, tdim, iseed,
                                   iref, algnd,
                                   &this->poi2bak(ipoin,tdim-1));

  if(this->param->iverb >= 3){
    printf("  -- END interpMetBack ipoin = %d ierro %d met = ",ipoin,ierro);
    dblAr1((this->idim*(this->idim+1))/2, this->met[ipoin]).print();
  }


  if(ierro != 0) return ierro;

  // No need to update back mesh if analytical metric. 
  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value) return 0;

  // Update other poi2baks.
  if(pdim == 1){
    int iedgb = this->poi2bak(ipoin,1-1);
    METRIS_ASSERT(iedgb >= 0 && iedgb < bak->nedge);
    int ifacb = bak->edg2fac[iedgb];
    METRIS_ASSERT(ifacb >= 0 && ifacb < bak->nface);
    this->poi2bak(ipoin,2-1) = ifacb;
  }
  if(this->get_tdim() > 2){
    int ifacb = this->poi2bak(ipoin,2-1);
    METRIS_ASSERT(ifacb >= 0 && ifacb < bak->nface);
    int itetb = bak->fac2tet[ifacb][0];
    METRIS_ASSERT(itetb >= 0 && itetb < bak->nelem);
    this->poi2bak(ipoin,3-1) = itetb;
  }

  return 0;
}

template<class MFT>
int Mesh<MFT>::interpMetBack0(int ipoi0,
                              int tdim, int iseed, 
                              int iref, 
                              const double*__restrict__ algnd,
                              int*__restrict__ ieleb){

  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
    this->met.getMetPhys(AsDeg::Pk,DifVar::None,this->met.getSpace(),
                         NULL,this->coord[ipoi0],this->met[ipoi0],NULL); 
  }else{

    const int iverb = this->param->iverb;

    if(this->idim == 3) METRIS_THROW_MSG(TODOExcept(), 
                                         "Metric interpolation in surface case")

    double barb[4];
    double coopr[3];
    bool ifnd = false;
    METRIS_ENFORCE(this->bak->idim == 2); // Otherwise put coop in a 3-sized 
    const int nnode = tdim == 1 ? edgnpps[this->curdeg] 
                    : tdim == 2 ? facnpps[this->curdeg] : tetnpps[this->curdeg];

    const intAr2& ent2poi = this->ent2poi(tdim);

    int pdim0 = this->get_tdim();
    int ibpo0 = this->poi2bpo[ipoi0];
    if(ibpo0 >= 0) pdim0 = this->bpo2ibi(ibpo0,1);
    const double *uvsrf = ibpo0 < 0 ? NULL : this->bpo2rbi[ibpo0];

    // First time around, skip ref mismatches
    // Second time, try localizing with no ref expectation. 
    for(int itry_ref = 0; itry_ref <= 1; itry_ref++){

      bool iskipped_lowdim = false;

      // Do two passes: in the first, we skip lower dim points 
      for(int iskiplow = 0; iskiplow <= 1; iskiplow++){
        if(iskiplow == 1 && !iskipped_lowdim) break;

        for(int ii = 0; ii < nnode; ii++){
          //int ipoin = this->fac2poi(iseed,ii);
          int ipoin = ent2poi(iseed,ii);
          int pdim  = this->getpoitdim(ipoin);
          if( (pdim > tdim) || (pdim < tdim && !iskipped_lowdim) ){
            if(iverb >= 3)
              printf(" - interpMetBack skip seed pt %d dim = %d > %d\n",
                     ipoin, pdim, tdim);
            if(pdim < tdim) iskipped_lowdim = true;
            continue;
          }

          *ieleb = poi2bak(ipoin, tdim-1);
          METRIS_ASSERT_MSG(*ieleb >= 0 && *ieleb < bak->nentt(tdim),"with tdim = "<<tdim);

          if(iverb >= 3){
            printf(" dump ipoin %d poi2bak = ",ipoin);
            for(int ii = 1; ii <= this->get_tdim(); ii++)
              printf(" %d : %d, ",ii,poi2bak(ipoin,ii-1));
            printf("\n");
          }

          int iref1 = iref;
          // If the point has less than tdim topo dim, we need to find a back seed
          // that has ref iref
          if(pdim < tdim){
            const intAr1& bakref = bak->ent2ref(tdim);
            if(iref >= 0 && iref != bakref[*ieleb]){
              // If first try, just skip this seed. 
              if(itry_ref == 0){
                if(iverb >= 3) printf (" - first try skip bad seed\n");
                continue;
              }else{
                // If second try, do localization without iref expectation. 
                // Do a posteriori check. 
                iref1 = -1;
              }
            }
          }

          bool dbgdist = false;

          if(iverb >= 3) 
            printf(" try locMesh with iref1 = %d iref %d ieleb %d pdim %d tdim %d"
              " bakref %d\n",iref1,iref,*ieleb,pdim,tdim,bak->ent2ref(tdim)[*ieleb]);

          CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bdeg == this->bak->curdeg){
            int ierro;

            #ifndef NDEBUG
            try{
            #endif

            if(this->idim == 2){
              if(tdim == 1){
                ierro = locMesh<2,1,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                          pdim0, uvsrf, iref1, algnd, coopr, barb);
              }else if(tdim == 2){
                ierro = locMesh<2,2,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                          pdim0, uvsrf, iref1, NULL , coopr, barb);
              }
            }else{
              if(tdim == 1){
                ierro = locMesh<3,1,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                          pdim0, uvsrf, iref1, algnd, coopr, barb);
              }else if(tdim == 2){
                ierro = locMesh<3,2,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                          pdim0, uvsrf, iref1, algnd, coopr, barb);
              }else{
                ierro = locMesh<3,3,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                          pdim0, uvsrf, iref1, NULL , coopr, barb);
              }
            }

            #ifndef NDEBUG
              double dist;
              if(this->idim == 2){
                dist = geterrl2<2>(this->coord[ipoi0],coopr);
              }else{
                dist = geterrl2<2>(this->coord[ipoi0],coopr);
              }
              if(sqrt(dist) >= 1){
                printf("## LARGE DISTANCE ? %f \n", dist);
                dbgdist = true;
                if(iverb >= 3){
                  int ipdbg = this->newpoitopo(-1,-1);
                  int ibdbg = this->template newbpotopo<0>(ipdbg,ipdbg);
                  for(int ii = 0; ii < this->idim; ii++) 
                    this->coord[ipdbg][ii] = this->coord(ipoi0,ii);
                  writeBackLinks("debug_poi2bak", *this);
                  writeMesh("debug_poi2bak", *this);
                  this->bpo2ibi(ibdbg,0) = -1;
                  this->set_nbpoi(this->nbpoi - 1);
                  this->set_npoin(this->npoin - 1);

                  if(dist > 1000){
                    printf("## EXIT HERE\n");
                    exit(1);
                  }
                }

              }
            #endif


            #ifndef NDEBUG
            }catch(const MetrisExcept &e){
              printf("## EXCEPTION THROWN IN LOCMESH, RERUN WITH PRINTS:\n");
              this->param->iverb = 10;

              int ipdbg = this->bak->newpoitopo(-1,-1);
              int ibdbg = this->bak->template newbpotopo<0>(ipdbg,ipdbg);
              for(int ii = 0; ii < this->idim; ii++) 
                this->bak->coord(ipdbg,ii) = this->coord(ipoi0,ii);

              printf("Try to localize coop %d = ",ipdbg);
              dblAr1(this->idim,this->coord[ipoi0]).print();
              dblAr1(this->idim,this->bak->coord[ipdbg]).print();
              writeMesh("debug-localization.meshb", *(this->bak));
              this->bak->bpo2ibi(ibdbg,0)  = -1;
              this->bak->set_npoin(this->bak->npoin-1);
              this->bak->set_nbpoi(this->bak->nbpoi-1);


              printf("WAIT HERE before throw\n");
              wait();
              if(this->idim == 2){
                if(tdim == 1){
                  ierro = locMesh<2,1,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                            pdim0, uvsrf, iref1, algnd, coopr, barb);
                }else if(tdim == 2){
                  ierro = locMesh<2,2,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                            pdim0, uvsrf, iref1, NULL, coopr, barb);
                }
              }else{
                if(tdim == 1){
                  ierro = locMesh<3,1,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                            pdim0, uvsrf, iref1, algnd, coopr, barb);
                }else if(tdim == 2){
                  ierro = locMesh<3,2,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                            pdim0, uvsrf, iref1, algnd, coopr, barb);
                }else{
                  ierro = locMesh<3,3,bdeg>(*(this->bak), ieleb, this->coord[ipoi0],
                                            pdim0, uvsrf, iref1, NULL, coopr, barb);
                }
              }

              printf("WAIT HERE before throw\n");
              wait();
              throw(e);
            }
            #endif

            // Check if projected is close in the metric, then keep it. 
            // Error LOC_ERR_ALLPOS can be caused by point actually outside domain
            // And conversely, loc might not error out but default to 
            // projection on boundary, could be an error if point too far
            // (e.g. loc stuck opposite side of hole)
            METRIS_ASSERT(*ieleb >= 0 && *ieleb < bak->nentt(tdim));
            METRIS_ASSERT((iref == bak->ent2ref(tdim)[*ieleb]) || iref == -1);

            
            const intAr2& bak2poi = bak->ent2poi(tdim);
            this->bak->met.getMetBary(AsDeg::Pk,
                                      DifVar::None,
                                      this->met.getSpace(),
                                      bak2poi[*ieleb],tdim,
                                      barb,this->met[ipoi0],NULL);

            double tang[3];
            double len;
            for(int ii = 0; ii < this->idim; ii++)
              tang[ii] = this->coord(ipoi0,ii) - coopr[ii];

            if(this->idim == 2){
              if(this->met.getSpace() == MetSpace::Log){
                len = getlenedg_log<2>(tang,this->met[ipoi0],100,1.0e-6);
              }else{
                len = getlenedg<2>(tang,this->met[ipoi0]);
              }
            }else{
              if(this->met.getSpace() == MetSpace::Log){
                len = getlenedg_log<3>(tang,this->met[ipoi0],100,1.0e-6);
              }else{
                len = getlenedg<3>(tang,this->met[ipoi0]);
              }
            }

            bool dowait = false;
            if(iverb >= 2 && ierro != 0 || iverb >= 3)
              printf("  - localization outside? len = %15.7e\n",len);
            if(iverb >= 2 && ierro != 0 || iverb >= 3){
              int nnmet = (this->idim*(this->idim+1))/2;
              printf(" using metl = ");
              dblAr1(nnmet,this->met[ipoi0]).print();
              printf(" in iele = %d bary ",*ieleb);
              dblAr1(tdim+1,barb).print();
              if(abs(barb[0]) + abs(barb[1]) + abs(barb[2]) > 100){
                printf("## VERY LARGE BARYCENTRIC COORDINATES?\n");
                dowait = true;
              }
            }
            if(len < 0.5){
              if(iverb >= 2 && ierro != 0 || iverb >= 3){
                printf("  -> len %15.7e < 0.5 keep w met = ",len);
                dblAr1( (this->idim*(this->idim+1))/2,this->met[ipoi0]).print();
              } 
              if(dbgdist && iverb >= 3){
                printf("Waiting check debug_poi2bak -> this point was kept despite large distance\n");
                //wait();
              }
              ierro = 0;
            }
            if(iverb >= 4 || len >= 0.5 && iverb >= 3){
              int ipdbg = this->bak->newpoitopo(-1,-1);
              int ibdbg = this->bak->template newbpotopo<0>(ipdbg,ipdbg);
              for(int ii = 0; ii < this->idim; ii++) 
                this->bak->coord[ipdbg][ii] = this->coord(ipoi0,ii);
              writeMesh("debug-localization.meshb", *(this->bak));
              this->bak->bpo2ibi(ibdbg,0) = -1;
              this->bak->set_npoin(this->bak->npoin - 1);
              this->bak->set_nbpoi(this->bak->nbpoi - 1);
              if(len >= 0.5) dowait = true;
            }

            if(dowait){
              printf("## WAIT HERE DUE TO POSSIBLE ERROR\n");
              int nnodb = this->bak->nnode(tdim);
              const intAr2& ent2pob = this->bak->ent2poi(tdim);
              printf("Localized in %d : ",*ieleb);
              intAr1(nnodb,ent2pob[*ieleb]).print();
              this->bak->setBasis(FEBasis::Lagrange);
              for(int ii = 0; ii < nnodb; ii++){
                int ipoib = ent2pob(*ieleb,ii);
                printf("%d : ",ipoib);
                for(int jj = 0; jj < this->idim; jj++) printf(" %23.15e ",
                  this->bak->coord(ipoib,jj));
                printf("\n");
              }
              printf("To loc: ");
              for(int jj = 0; jj < this->idim; jj++) printf(" %23.15e ",
                                                        this->coord(ipoi0,jj));
              printf("In basis ");
              if(this->bak->getBasis() == FEBasis::Lagrange) printf(" Lagrange \n");
              else printf(" Bézier \n");


              bool iinva;
              double ccoef[tetnpps[this->idim*(bdeg-1)]];
              if(this->idim == 2){
                getsclccoef<2,2,bdeg>(*(this->bak),*ieleb,NULL,ccoef,&iinva);
              }else{
                getsclccoef<3,3,bdeg>(*(this->bak),*ieleb,NULL,ccoef,&iinva);
              }
              printf("Element invalid ? %d \n",iinva);
              wait();
            }

            if(ierro == 0){
              ifnd = true;
              if(iverb >= 4){
                printf("   - interpMetBack loc in %d dim %d bary = ",*ieleb,tdim);
                dblAr1(tdim+1,barb).print();
                printf("   - metric = ");
                dblAr1( (this->idim*(this->idim+1))/2,this->met[ipoi0]).print();
              }
            }

          }}CT_FOR1(bdeg);
          if(ifnd) return 0;
        }// for ii 
      } // for iskiplow
    }// for itry_ref
    if(!ifnd) return 1;
  }

  // Point has been found 
  // Update all the poi2baks of higher topological dimension 



  return 0;
}

template class Mesh<MetricFieldFE>;
template class Mesh<MetricFieldAnalytical>;


} // End namespace
