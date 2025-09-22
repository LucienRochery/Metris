//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "Mesh.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../metris_constants.hxx"
#include "../ho_constants.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_topo.hxx"
#include "../low_geo/lenedg.hxx"
#include "../low_geo/misc.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/ccoef.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../utils/aux_misc.hxx"

#include "../Localization/msh_localization.hxx"
#include "../Localization/msh_interpFrontBack.hxx"


namespace Metris{


template<class MFT>
int Mesh<MFT>::interpMetBack(int ipoin){
  GETVDEPTH(this->param);
  METRIS_ASSERT(ipoin >= 0 && ipoin < this->npoin);
  
  int tdim = this->getpoitdim(ipoin);
  METRIS_ASSERT_MSG(tdim > 0, "Interpolate corners manually");
  
  int iseed = this->poi2ent(ipoin,0);
  // Control point case
  if(iseed < -1) iseed = - iseed - 2;
  METRIS_ASSERT(this->poi2ent(ipoin,1) == tdim);
  METRIS_ASSERT(iseed >= 0 && iseed < this->nentt(tdim));
  METRIS_ASSERT(!isdeadent(iseed,this->ent2poi(tdim)));
  
  int iref = this->ent2ref(tdim)[iseed];
  METRIS_ASSERT(iref >= 0);

  CPRINTF2("-- START interpMetBack ipoin {} pdim {} iseed {} iref {} poi2bak {}\n",
           ipoin,tdim,iseed,iref,this->poi2bak[ipoin]);

  // Get algnd
  double algnd[3];
  this->get_algnd(ipoin, algnd);

  // First try seeding point with itself, if already initialized. 
  if(this->poi2bak[ipoin] >= 0){
    if(this->interpMetBack00(ipoin, tdim, iref, ipoin, algnd) == 0) return 0;
    CPRINTF1(" - self-seeded interpmet back failed, retry with seed {} \n", iseed);
  }

  return this->interpMetBack(ipoin,tdim,iseed,iref,algnd);
}

template int Mesh<MetricFieldFE>::interpMetBack(int ipoin);
template int Mesh<MetricFieldAnalytical>::interpMetBack(int ipoin);



// tdim is the dimension of iseed (front elt), AND of the point !
// iref and algnd only necessary if doing boundary localization. 
template<class MFT>
int Mesh<MFT>::interpMetBack(int ipoin, int tdim, int iseed, 
                             int iref, const double* algnd){
  METRIS_ASSERT_MSG(tdim == this->getpoitdim(ipoin) || this->getpoitdim(ipoin) == 0,
    "seed is dim {} point is {} ipoin = {}", tdim, this->getpoitdim(ipoin), ipoin);
  METRIS_ASSERT_MSG(iseed >= 0 && iseed < this->nentt(tdim), "interpMetBack provided invalid seed index {}", iseed);
  METRIS_ASSERT_MSG(!isdeadent(iseed,tdim == 1 ? this->edg2poi :
                                     tdim == 2 ? this->fac2poi : this->tet2poi),
                    "Dead seed passed to interpMetBack");

  GETVDEPTH(this->param);

  int pdim = this->getpoitdim(ipoin);
  METRIS_ASSERT(pdim == tdim);

  CPRINTF1("-- START interpMetBack ipoin = {} iseed = {} tdim {} \n",ipoin,
           iseed,tdim);

  METRIS_ASSERT_MSG(tdim == this->idim || (algnd != NULL && iref >= 0),
    "Boundary dim {} but algnd = {} and iref = {}", tdim, (void*) algnd, iref);

  METRIS_ASSERT_MSG(ipoin >= 0 && ipoin < this->npoin, 
    "interpMetBack ipoin out of bounds {} < ? {}", ipoin, this->npoin);

  int ierro = this->interpMetBack0(ipoin, tdim, iseed, iref, algnd);

  if(DOPRINTS1()){
    CPRINTF1("-- END interpMetBack ipoin = {} ierro {} met = {}\n",
             ipoin,ierro,dblAr1((this->idim*(this->idim+1))/2, this->met[ipoin]));
  }

  return ierro;
}

template int Mesh<MetricFieldFE>::interpMetBack(int ipoin, int tdim, int iseed, 
                                                 int iref, const double* algnd);
template int Mesh<MetricFieldAnalytical>::interpMetBack(int ipoin, int tdim, int iseed, 
                                                 int iref, const double* algnd);

template<class MFT>
int Mesh<MFT>::interpMetBack0(int ipoi0,int tdim, int iseed, int iref, 
                              const double*__restrict__ algnd){

  GETVDEPTH(this->param)

  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
    this->met.getMetPhys(DifVar::None,this->met.getSpace(),
                         this->coord[ipoi0],this->met[ipoi0],NULL); 
    return 0;
  }else{

    //if(this->idim == 3) METRIS_THROW_MSG( 
    //                                     "Metric interpolation in surface case")

    const int nnode = tdim == 1 ? getnnod1(this->curdeg) 
                    : tdim == 2 ? getnnod2(this->curdeg) : getnnod3(this->curdeg);

    const intAr2& ent2poi = this->ent2poi(tdim);
    METRIS_ASSERT_MSG(!isdeadent(iseed,ent2poi),"Dead seed passed to interpMetBack");

    // First time around, skip ref mismatches
    // Second time, try localizing with no ref expectation. 
    for(int itry_ref = 0; itry_ref <= 1; itry_ref++){

      bool iskipped_lowdim = false;

      // Do two passes: in the first, we skip lower dim points 
      for(int iskiplow = 0; iskiplow <= 1; iskiplow++){
        if(iskiplow == 1 && !iskipped_lowdim) break;

        CPRINTF2(" - interpMetBack pass {}/{} (reject 0/accept 1: ref mismatch/low dim)\n",itry_ref,iskiplow);

        for(int ii = 0; ii < nnode; ii++){
          //int ipoin = this->fac2poi(iseed,ii);
          int ipseed = ent2poi(iseed,ii);
          int pdim  = this->getpoitdim(ipseed);
          //if( (pdim > tdim) || (pdim < tdim && !iskipped_lowdim) ){
          if( (pdim > tdim) || (pdim < tdim && iskiplow == 0) ){
            CPRINTF2(" - skip seed pt {} dim = {} != {}\n",
                     ipseed, pdim, tdim);
            if(pdim < tdim) iskipped_lowdim = true;
            continue;
          }

          if(this->interpMetBack00(ipoi0, tdim, iref, ipseed, algnd) == 0) return 0;

        }// for ii 
      } // for iskiplow
    }// for itry_ref
    return 1;
  }

  // Point has been found 
  // Update all the poi2baks of higher topological dimension 



  return 0;
}

template int Mesh<MetricFieldFE>::interpMetBack0(int ipoi0,int tdim, int iseed, 
                                     int iref, const double*__restrict__ algnd);
template int Mesh<MetricFieldAnalytical>::interpMetBack0(int ipoi0,int tdim, 
                          int iseed, int iref, const double*__restrict__ algnd);


// ipoi0 to locate of dimension pdim0 with ref constraint iref (or < 0)
// ipseed in front mesh has been interpMetBack's in the past
// algnd is tangent or normal for localization, as usual. 
template<class MFT>
int Mesh<MFT>::interpMetBack00(int ipoi0, int tdim, int iref, int ipseed,
                               const double*__restrict__ algnd){
  GETVDEPTH(this->param);

  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
    this->met.getMetPhys(DifVar::None,this->met.getSpace(),
                         this->coord[ipoi0],this->met[ipoi0],NULL); 
    return 0;
  }

  CPRINTF1("-- START interpMetBack00 ipoi0 {} tdim {} iref {} ipseed {}\n",
           ipoi0,tdim,iref,ipseed);

  int ieleb = poi2bak[ipseed];
  if(ieleb < 0){
    CPRINTF2(" - skip seed {} < 0\n",ieleb);
    return 2; // Happens when called from the cavity operator.
  }

  CPRINTF2(" - init ieleb = {}\n",ieleb);

  int pdim_seed;

  // A point can be self-seeded; in that case, the ieleb is already correct dim
  if(ipseed == ipoi0) goto ieleb_initialized;

  pdim_seed = this->getpoitdim(ipseed);

  if(pdim_seed > tdim){
    CPRINTF2(" - seed has dim {} > {} = point dim, skip\n",pdim_seed,tdim);
    return 3;
  }

  CPRINTF1(" - seed point dim {} initial bak seed {}\n",pdim_seed,ieleb);

  if(pdim_seed == tdim){
    // Easiest case, ieleb is already the one, and we just need to check ref.
    if(iref >= 0 && this->bak->ent2ref(tdim)[ieleb] != iref){
      CPRINTF1(" # provided seed has ref {} != {} \n",this->bak->ent2ref(tdim)[ieleb],iref);
      return 4;
    }
    goto ieleb_initialized;
  }

  // Otherwise, the seed is dim < tdim: it can have several refs elts of dim tdim.
  if(pdim_seed == 0){
    // This is the easiest case, because we get a link to a vertex (corner)
    // in the back mesh. 
    int ibpob = this->bak->poi2ebp(ieleb, tdim, -1, iref);
    if(ibpob < 0){
      CPRINTF1(" # could not find seed element dim {} ref {} from point {}\n",
               tdim, iref, poi2bak[ipseed]);
      return 5;
    }
    ieleb = this->bak->bpo2ibi(ibpob,2);
    METRIS_ASSERT(this->bak->ent2ref(tdim)[ieleb] == iref);
    CPRINTF1(" - corner seed tdim {} element {}\n",tdim, ieleb);
    goto ieleb_initialized;
  }
  // If the seed is ref unconstrained, we will simply use edg2fac and fac2tet.
  // This is either if no ref constraint, or if single ref of target dim.
  if(iref < 0 || (tdim == 1 && this->CAD.ncaded == 1)
              || (tdim == 2 && this->CAD.ncadfa == 1)
              || (tdim == 3 && this->ndomn      == 1)){
    if(pdim_seed == 1){
      ieleb = this->bak->edg2fac[ieleb];
      CPRINTF1(" - seed edge -> face = {}\n",ieleb);
    }
    if(pdim_seed <= 2 && tdim == 3){
      ieleb = this->bak->fac2tet(ieleb,0);
      if(ieleb < 0) ieleb = this->bak->fac2tet(ieleb,1);
      CPRINTF1(" - seed face -> tetra = {}\n",ieleb);
    }
    METRIS_ASSERT(ieleb >= 0);
    goto ieleb_initialized;
  }

  // Last case is pdim_seed >= 1, iref >= 0. If tdim is boundary, we have 
  // poi2ebp. We can use this on a vertex of the back seed. 
  if(this->isboundary_tdim(tdim)){
    int nnode = getnnode(pdim_seed, this->bak->curdeg);
    for(int inode = 0; inode < nnode; inode++){
      int ipoib = this->bak->ent2poi(pdim_seed)(ieleb, inode);
      int ibpob = this->bak->poi2ebp(ipoib, tdim, -1, iref);
      if(ibpob < 0) continue;
      CPRINTF1(" - back seed elt {} vertex {} = {} points to back elt {} of ref {}\n",
               ieleb, inode, ipoib, this->bak->bpo2ibi(ibpob, 2), iref);
      ieleb = this->bak->bpo2ibi(ibpob, 2);
      METRIS_ASSERT(ieleb >= 0 && ieleb < this->bak->nentt(this->get_tdim()));
      goto ieleb_initialized;
    }
    return 6;
  }else{
    static int nwarnprt = 0;
    if(nwarnprt++ < 10) printf("## IMPLEMENT BACK SEED INIT IN MULTI DOMAIN MESH");
    return 7;
  }

ieleb_initialized:

  METRIS_ASSERT_MSG(ieleb >= 0 && ieleb < bak->nentt(tdim),
    "with tdim = {} got ieleb = {} ipseed = {}", tdim, ieleb, ipseed);

  CPRINTF2(" - using ipseed {} bak elt seed = {} dim {}\n", ipseed, ieleb, tdim);


  int ierro;
  double barb[4], coopr[3];

  double *uvsrf = NULL;
  if(tdim < this->get_tdim()){
    int ibpo0 = this->poi2bpo[ipoi0];
    METRIS_ASSERT(ibpo0 >= 0);
    if(this->bpo2ibi(ibpo0,1) == 0) ibpo0 = this->bpo2ibi(ibpo0,3);
    METRIS_ASSERT(ibpo0 >= 0);
    uvsrf = this->bpo2rbi[ibpo0];
  }

  #ifndef NDEBUG
  try{
  #endif
  CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bdeg == this->bak->curdeg){
    INCVDEPTH(this->param);
    if(this->idim == 2){
      METRIS_ASSERT(tdim <= 2);
      if(tdim == 1){
        ierro = locMesh<2,1,bdeg>(*(this->bak), &ieleb, this->coord[ipoi0],
                                  tdim, uvsrf, iref, algnd, coopr, barb);
      }else if(tdim == 2){
        ierro = locMesh<2,2,bdeg>(*(this->bak), &ieleb, this->coord[ipoi0],
                                  tdim, uvsrf, iref, NULL , coopr, barb);
      }
    }else{
      if(tdim == 1){
        ierro = locMesh<3,1,bdeg>(*(this->bak), &ieleb, this->coord[ipoi0],
                                  tdim, uvsrf, iref, algnd, coopr, barb);
      }else if(tdim == 2){
        ierro = locMesh<3,2,bdeg>(*(this->bak), &ieleb, this->coord[ipoi0],
                                  tdim, uvsrf, iref, algnd, coopr, barb);
      }else{
        ierro = locMesh<3,3,bdeg>(*(this->bak), &ieleb, this->coord[ipoi0],
                                  tdim, uvsrf, iref, NULL , coopr, barb);
      }
    }
  }}CT_FOR1(bdeg);


  CPRINTF1(" - locMesh return ierro {} ieleb = {} tdim {} \n",ierro, ieleb, tdim);

  #ifndef NDEBUG
  }catch(const MetrisExcept &e){
    PRINTF("## EXCEPTION THROWN IN LOCMESH, RERUN WITH PRINTS:\n");
    this->param->iverb   = 10;
    this->param->ivdepth = 20;

    int ipdbg = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    int ibdbg = this->bak->newbpotopo(Vertex{ipdbg},0,ipdbg);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdbg,ii) = this->coord(ipoi0,ii);

    int ipdb2 = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    this->bak->newbpotopo(Vertex{ipdb2},0,ipdb2);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdb2,ii) = coopr[ii];


    PRINTF("Try to localize coop {} = {}, bak ipdbg = {}\n",ipdbg,dblAr1(this->idim,this->coord[ipoi0]),dblAr1(this->idim,this->bak->coord[ipdbg]));
    writeMesh("debug-localization.meshb", *(this->bak));
    this->bak->bpo2ibi(ibdbg,0)  = -1;
    this->bak->killpoint(ipdbg);
    this->bak->killpoint(ipdb2);

    writeMesh("debug-front.meshb", *(this));

    PRINTF("WAIT HERE before throw\n");
    wait();
    throw(e);
  }// catch
  #endif

  if(ierro != 0 && DOPRINTS2() && this->param->dbgfull){
    int ipdbg = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    int ibdbg = this->bak->newbpotopo(Vertex{ipdbg},0,ipdbg);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdbg,ii) = this->coord(ipoi0,ii);

    int ipdb2 = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    this->bak->newbpotopo(Vertex{ipdb2},0,ipdb2);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdb2,ii) = coopr[ii];


    PRINTF("Try to localize coop {} = {}, bak ipdbg = {} = {}\n",ipdbg,dblAr1(this->idim,this->coord[ipoi0]), ipdbg, dblAr1(this->idim,this->bak->coord[ipdbg]));
    writeMesh("debug-localization.meshb", *(this->bak));
    this->bak->bpo2ibi(ibdbg,0)  = -1;
    this->bak->killpoint(ipdbg);
    this->bak->killpoint(ipdb2);
    writeMesh("debug-front.meshb", *(this));
  }

  // Check if projected is close in the metric, then keep it. 
  // Error LOC_ERR_ALLPOS can be caused by point actually outside domain
  // And conversely, loc might not error out but default to 
  // projection on boundary, could be an error if point too far
  // (e.g. loc stuck opposite side of hole)
  METRIS_ASSERT(ieleb >= 0 && ieleb < bak->nentt(tdim));
  METRIS_ASSERT_MSG((iref == bak->ent2ref(tdim)[ieleb]) || iref == -1,
    "got iref = {} ieleb = {} tdim {} ieleb ref {}", iref, ieleb, tdim, bak->ent2ref(tdim)[ieleb]);


  
  const intAr2& bak2poi = bak->ent2poi(tdim);
  this->bak->met.getMetBary(AsDeg::Pk,
                            DifVar::None,
                            this->met.getSpace(),
                            bak2poi[ieleb],tdim,
                            barb,this->met[ipoi0],NULL);

  // In release mode, don't check the length in case the localization did not
  // return error. getlenedg is expensive.
  #ifdef NDEBUG
  if(ierro == 0){
    if(DOPRINTS2()){
      CPRINTF2(" - interpMetBack loc in {} dim {} bary = ",ieleb,tdim);
      dblAr1(tdim+1,barb).print();
      CPRINTF2(" - metric = ");
      dblAr1( (this->idim*(this->idim+1))/2,this->met[ipoi0]).print();
    }
    this->poi2bak[ipoi0] = ieleb;
    return 0;
  }
  #endif

  // Only in debug mode:
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

  if(DOPRINTS2()){
    CPRINTF2("- localization outside? len = {:15.7e} w tang {}\n",len,dblAr1(this->idim,tang));
    int nnmet = (this->idim*(this->idim+1))/2;
    CPRINTF2(" using metl = {}\n",dblAr1(nnmet,this->met[ipoi0]));
    CPRINTF2(" in iele = {} bary {}\n",ieleb,dblAr1(tdim+1,barb));
    double sum = abs(barb[0]) + abs(barb[1]);
    if(tdim >= 2) sum += abs(barb[2]);
    if(tdim >= 3) sum += abs(barb[3]);
    if(sum > 100){
      CPRINTF2("## VERY LARGE BARYCENTRIC COORDINATES?\n");
    }
    intAr2 &ent2pob = this->bak->ent2poi(tdim);
    for(int ii = 0; ii < tdim + 1; ii++){
      CPRINTF2("vertex {} metric = {}\n",ent2pob(ieleb,ii),
               dblAr1(nnmet,this->bak->met[ent2pob(ieleb,ii)]));
    }
  }
  if(len < 0.5){
    CPRINTF2("-> len {:15.7e} < 0.5 keep w met = {}\n",len,dblAr1((this->idim*(this->idim+1))/2,this->met[ipoi0]));
    ierro = 0;
  }else{
    CPRINTF1("# Large length {} \n",len);
    CPRINTF2("- proj point {}\n",dblAr1(this->idim,coopr));
    CPRINTF2("- loc  point {}\n",dblAr1(this->idim, this->coord[ipoi0]));
    ierro = 1;
  }

  if(DOPRINTS2() && len >= 0.5){
    int ipdbg = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    int ibdbg = this->bak->newbpotopo(Vertex{ipdbg},0,ipdbg);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdbg,ii) = this->coord(ipoi0,ii);

    int ipdb2 = this->bak->newpoitopo(PointType::Vertex,-1,-1);
    this->bak->newbpotopo(Vertex{ipdb2},0,ipdb2);
    for(int ii = 0; ii < this->idim; ii++) 
      this->bak->coord(ipdb2,ii) = coopr[ii];

    const int nnmet = (this->idim*(this->idim + 1))/2;
    double metl[6];

    for(int ii = 0; ii < nnmet; ii++)
      metl[ii] = this->met(ipoi0,ii);

    if(this->met.getSpace() == MetSpace::Exp){
      if(this->idim == 2){
        getlogmet_inp<2>(metl);
      }else{
        getlogmet_inp<3>(metl);
      }
    }

    for(int ii = 0; ii < nnmet; ii++){
      this->bak->met(ipdbg,ii) = metl[ii];
      this->bak->met(ipdb2,ii) = metl[ii];
    }

    writeMesh("debug-localization.meshb", *(this->bak));
    this->bak->met.writeMetricFile("debug-localization.solb");
    this->bak->bpo2ibi(ibdbg,0) = -1;
    this->bak->killpoint(ipdbg);
    this->bak->killpoint(ipdb2);
    writeMesh("debug-front.meshb", *(this));
  }
  if(ierro == 0){
    this->poi2bak[ipoi0] = ieleb;
  }

  return ierro > 0 ? 10+ierro : 0;
}
template int Mesh<MetricFieldFE>::interpMetBack00(int ipoi0, int tdim, 
                          int iref, int ipseed,const double*__restrict__ algnd);
template int Mesh<MetricFieldAnalytical>::interpMetBack00(int ipoi0, int tdim, 
                          int iref, int ipseed,const double*__restrict__ algnd);

} // End namespace
