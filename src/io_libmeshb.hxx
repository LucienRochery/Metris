//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __IO_LIBMESHB__
#define __IO_LIBMESHB__


#include <string>
#include "../libs/libmeshb.hxx"
#include "Mesh/MeshFwd.hxx"
#include "types.hxx"


namespace Metris{

class MshCavity;
// For writeField 
enum class SolTyp{P0Elt,CG};

void test_read();

std::string correctExtension_meshb(const std::string &s);
std::string correctExtension_solb(const std::string &s);
std::string correctExtension_egads(const std::string &s);

// Wrapper handling exceptions
// For writing
template<int rwtyp>
int64_t MetrisOpenMeshFile(std::string name, int meshDim);
// For reading
template<int rwtyp>
int64_t MetrisOpenMeshFile(std::string name, int *meshDim);

//void readMeshFile(int64_t libIdx, MeshBase &msh, int iprt);

//int iniMesh(cxxopts::ParseResult &paropt, int usrMinDeg, int usrMaxDeg, Mesh &msh, Mesh &bak);
//void iniMesh(MetrisOptions &opt, int usrMinDeg, int usrMaxDeg, Mesh &msh);

void writeMeshCavity(std::string meshName, MeshBase &msh, const MshCavity& cav, int ithread = 0);
void writeMesh(std::string meshName, MeshBase &msh, bool ivolonly = false,
               int nedg0 = 0, int nfac0 = 0, int nele0 = 0);

template<class MFT>
void writeBackLinks(std::string solName, Mesh<MFT>& msh);


void debugInveval(std::string meshName_, MeshBase &msh, int tdim, int* ent2pol, double *coop);

#if 0
void writeMesh(std::string meshName, int ideg, int ilag,
               int npoin, dblAr2 &coord,
               int nelem, intAr2 &tet2poi, intAr1 &tet2ref,
               int nface, intAr2 &fac2poi, intAr1 &fac2ref,
               int nedge, intAr2 &edg2poi, intAr1 &edg2ref);
void writeMeshVecs(std::string meshName, MeshBase &msh, const dblAr2 &poi2vec);
#endif

void writeField(std::string outname, const MeshBase &msh, SolTyp stype, dblAr1 &rfld, int ndim = 1);

void writeEdgesLengths(const MeshBase &msh, std::string outnroot, 
                       intAr2 &edg2poi, const dblAr1 &rlened);


//void writeMetric(std::string metName, const Mesh &msh);
//void writeExpMetric(std::string metName, Mesh &msh);

} // End namespace
#endif