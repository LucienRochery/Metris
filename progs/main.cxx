//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "main_adap.hxx"
#include "metris_options.hxx"
/*
(tet2poi|tet2tet)\((iele[m2]),([a-Z0-9]*)\)
(tet2poi|tet2tet)\((iele[m2]),(lno(fa|ed)3\[i(fa|ed)2?\]\[[0-9a-Z]\])\)
$1[$2][$3]
*/



int main(int argc, char** argv){



  int icod;

    //Mesh msh, bak;
    try{
      icod = Metris::main_metris(argc, argv);
      //icod = main_metris(argc, argv, msh, bak);
    }catch(const Metris::MetrisExcept &e){
      fmt::print("\n################################################################\n");
      fmt::print(stderr,"## MAIN_METRIS THROWS EXCEPTION:\n");
      fmt::print(stderr,"## Message: {}\n",e.message);

    #ifndef NO_BOOST_EXCEPT
      if(std::string const * ms=boost::get_error_info<excMessage>(e) )
        std::cout<<"## Message: "<<*ms;
      if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
        std::cerr << "## Call stack: \n" << *tr;
    #endif
    }

  return icod;
}


