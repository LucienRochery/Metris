
  int ipcor = 168;
  printf("## DEBUG print corner ibpois %d \n",ipcor);
  for(int ibpoi = msh.poi2bpo[ipcor]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
    printf(" ibpoi %d : ",ibpoi);
    intAr1(nibi, msh.bpo2ibi[ibpoi]).print();
    printf(" t/u,v = %e %e \n",msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
    int ientt = msh.bpo2ibi(ibpoi,2);
    int tdimn = msh.bpo2ibi(ibpoi,1);
    if(tdimn == 0) continue;
    int iref = msh.ent2ref(tdimn)[ientt];
    printf(" ientt %d dim %d ref %d\n", ientt,tdimn,iref);
  }


  intAr2 fac2eo0(msh.CAD.ncadfa, msh.CAD.ncaded);
  std::map<ego,int> cde2ref;
  std::map<ego,int> cdf2ref;

  for(int ii = 0; ii < msh.CAD.ncaded; ii++){
    ego edge = msh.CAD.cad2edg[ii];
    cde2ref[edge] = ii;
  }
  for(int ii = 0; ii < msh.CAD.ncadfa; ii++){
    ego face = msh.CAD.cad2fac[ii];
    cdf2ref[face] = ii;
  }

  printf("## DEBUG print out faces and edges orientations\n");
  for(int icadfa = 0; icadfa < msh.CAD.ncadfa; icadfa++){
    ego face = msh.CAD.cad2fac[icadfa], *lchild, geom; 
    int ireff = cdf2ref[face];
    int oclass,mtype,nchild,*senses;
    int ierro = EG_getTopology(face,&geom,&oclass,&mtype,NULL,
                               &nchild,&lchild,&senses);
    METRIS_ENFORCE(ierro == EGADS_SUCCESS);
    printf("face %d nchild = %d\n",icadfa,nchild);
    for(int ii = 0; ii < nchild; ii++){
      ego child = lchild[ii];
      printf("child mtype = %d oclass %d sens %d\n", child->mtype, child->oclass,
        senses[ii]);
    
      int oclass2,mtype2,nchild2,*senses2;
      ego *lchild2, geom2; 
      ierro = EG_getTopology(child,&geom2,&oclass2,&mtype2,NULL,
                             &nchild2,&lchild2,&senses2);
      for(int jj = 0; jj < nchild2; jj++){
        ego child2 = lchild2[jj];
        int irefe = cde2ref[child2];
        printf("  child mtype %d oclass %d ref %d sens %d\n", child2->mtype, child2->oclass,
          irefe,senses2[jj]);
        fac2eo0(ireff, irefe) = senses2[jj] * senses[ii];
      }
    }
  }

  for(int icaded = 0; icaded < msh.CAD.ncaded; icaded++){
    ego edge = msh.CAD.cad2edg[icaded];
    double range[2];
    int periodic;
    EG_getRange(edge, range, &periodic);
    printf(" edge ref %d range %e %e periodic %d\n",icaded, range[0], range[1],
      periodic);
  }

  intAr2 fac2eor(msh.CAD.ncadfa, msh.CAD.ncaded);
  fac2eor.fill(0);
  int nsucc = 0, nfail = 0;
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge, msh.edg2poi)) continue;

    int irefe = msh.edg2ref[iedge];

    //if(irefe != 11) continue;

    int ip1 = msh.edg2poi(iedge,0);
    int ip2 = msh.edg2poi(iedge,1);

    int ifac1 = msh.edg2fac[iedge];
    int ied1 = getedgfac(msh, ifac1, ip1, ip2);
    int ifac2 = msh.fac2fac(ifac1, ied1);
    int ied2 = getedgfac(msh, ifac2, ip1, ip2);

    int lface[2] = {ifac1, ifac2};
    int led[2] = {ied1, ied2};
    for(int ii = 0; ii < 2; ii++){
      int iface = lface[ii];
      int ied = led[ii];
      int iref = msh.fac2ref[iface];
      //if(iref != 8) continue;

      int ip1 = msh.fac2poi(iface,ied);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip3 = msh.fac2poi(iface,lnoed2[ied][1]);

      int ibf1 = msh.poi2ebp(ip1, 2, iface, iref);
      int ibf2 = msh.poi2ebp(ip2, 2, iface, iref);
      int ibf3 = msh.poi2ebp(ip3, 2, iface, iref);

      int ibe1 = msh.poi2ebp(ip2, 1, iedge, irefe);
      int ibe2 = msh.poi2ebp(ip3, 1, iedge, irefe);


      double t1 = msh.bpo2rbi(ibe1, 0);
      double t2 = msh.bpo2rbi(ibe2, 0);
      int isens = 1;
      if(t1 >= t2) isens = -1;

      //printf("iface %d ip %d %d %d = %d %d %d\n",iface,ip1,ip2,ip3,
      //  msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2));
      //printf("  ip %d ibf %d (u,v) = %e %e\n",ip1, ibf1, msh.bpo2rbi(ibf1,0), msh.bpo2rbi(ibf1,1));
      //printf("  ip %d ibf %d (u,v) = %e %e\n",ip2, ibf2, msh.bpo2rbi(ibf2,0), msh.bpo2rbi(ibf2,1));
      //printf("  ip %d ibf %d (u,v) = %e %e\n",ip3, ibf3, msh.bpo2rbi(ibf3,0), msh.bpo2rbi(ibf3,1));
      //printf("  ip %d ibe %d t = %e\n", ip2, ibe1, t1);
      //printf("  ip %d ibe %d t = %e\n", ip3, ibe2, t2);


      double meas = det2_vdif(msh.bpo2rbi[ibf2], msh.bpo2rbi[ibf1],
                              msh.bpo2rbi[ibf3], msh.bpo2rbi[ibf1]);
      meas *= isens;
      int iori = meas > 0 ? 1 : -1;
      //printf("  meas = %e isens %d ori = %d\n",meas/isens,isens,iori);
      if(fac2eor(iref, irefe) == 0){
        fac2eor(iref, irefe) = iori;
        printf(" face %5d ref %3d find ori %3d with edge %5d ref %3d EGADS: %3d\n",
          iface, iref, iori, iedge, irefe, fac2eo0(iref, irefe));
      }else if(fac2eor(iref,irefe) != iori){
        printf("## FACE %d ref %d find ori %d with edge %d ref %d, had %d\n",
          iface, iref, iori, iedge, irefe, fac2eor(iref,irefe));
        nfail++;
      }else{
        nsucc++;
      }

      bool doprints = false;

      int icode;
      double uvm[2], uvp[2];
      ego face = msh.CAD.cad2fac[iref];
      ego edge = msh.CAD.cad2edg[irefe];

      icode = EG_getEdgeUV(face, edge, 1, t1, uvp);
      icode = EG_getEdgeUV(face, edge,-1, t1, uvm);
      //printf("  t1 uvp %e %e uvm %e %e uv stored %e %e\n", 
      //       uvp[0], uvp[1], uvm[0], uvm[1], msh.bpo2rbi(ibf2,0), msh.bpo2rbi(ibf2,1));
      //if(geterrl2<2>(uvp, msh.bpo2rbi[ibf2]) > 1.0e-6){
      //  printf("## +1 Large uv diff:");
      //  doprints = true;
      //}
      //if(geterrl2<2>(uvm, msh.bpo2rbi[ibf2]) > 1.0e-6){
      //  printf("## -1 Large uv diff:");
      //  doprints = true;
      //}
      if(geterrl2<2>(uvm, uvp) > 1.0e-6){
        printf("## -+ Large uv diff:");
        doprints = true;
      }
      if(doprints){
        printf("  t1 uvp %e %e uvm %e %e uv stored %e %e\n", 
               uvp[0], uvp[1], uvm[0], uvm[1], msh.bpo2rbi(ibf2,0), msh.bpo2rbi(ibf2,1));
        printf("iface %d ip %d %d %d = %d %d %d\n",iface,ip1,ip2,ip3,
          msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2));
        printf("  ip %d ibf %d (u,v) = %e %e\n",ip1, ibf1, msh.bpo2rbi(ibf1,0), msh.bpo2rbi(ibf1,1));
        printf("  ip %d ibf %d (u,v) = %e %e\n",ip2, ibf2, msh.bpo2rbi(ibf2,0), msh.bpo2rbi(ibf2,1));
        printf("  ip %d ibf %d (u,v) = %e %e\n",ip3, ibf3, msh.bpo2rbi(ibf3,0), msh.bpo2rbi(ibf3,1));
        printf("  ip %d ibe %d t = %e\n", ip2, ibe1, t1);
        printf("  ip %d ibe %d t = %e\n", ip3, ibe2, t2);

        double result[18], resultp[18], resultm[18];
        icode = EG_evaluate(face, msh.bpo2rbi[ibf2], result);
        icode = EG_evaluate(face, uvp, resultp);
        icode = EG_evaluate(face, uvm, resultm);

        double err1 = geterrl2<3>(result, msh.coord[ip2]);
        double errp = geterrl2<3>(resultp, msh.coord[ip2]);
        double errm = geterrl2<3>(resultm, msh.coord[ip2]);

        printf(" # coord error from stored %e +1 %e -1 %e\n",err1,errp,errm);

      }
    }
  }
  printf("nsucc %d nfail %d\n",nsucc, nfail);
