//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_PRINT_MACROS__
#define __METRIS_PRINT_MACROS__

#include <string>
#include "aux_exceptions.hxx"
#include "aux_utils.hxx"
#include "MetrisRunner/MetrisParameters.hxx"

#define INCVDEPTH(has_param)  DepthCounter dc__(true);\
const char* spaces_string__ = dc__.getSpaces();\
const int iverb__ = has_param.param->iverb;\
const int ivdepth__ = has_param.param->ivdepth;


#define GETVDEPTH(has_param) DepthCounter dc__(false);\
const char* spaces_string__ = dc__.getSpaces();\
const int iverb__ = has_param.param->iverb;\
const int ivdepth__ = has_param.param->ivdepth;

//std::string spaces_string__ = dc__.getSpaces();\

// ##__VA_ARGS__ deletes trailing comma if __VA_ARGS__ empty. C++20 has more 
// elegant solutions but this should be most portable.
// Two verbosity levels.
#define DOPRINTS1() (dc__.getDepth() <= ivdepth__  && iverb__ >= 1)
#define DOPRINTS2() (dc__.getDepth() <= ivdepth__  && iverb__ >= 2)
#define MPRINTF(fmt,...) printf("%s" fmt, spaces_string__, ##__VA_ARGS__);
#define CPRINTF1(fmt,...) if(DOPRINTS1()){printf("%s" fmt, spaces_string__, ##__VA_ARGS__);}
#define CPRINTF2(fmt,...) if(DOPRINTS2()){printf("%s" fmt, spaces_string__, ##__VA_ARGS__);}
//#define CPRINTF1(fmt,...) if(DOPRINTS1()){printf("%s" fmt, spaces_string__.c_str(), ##__VA_ARGS__);}
//#define CPRINTF2(fmt,...) if(DOPRINTS2()){printf("%s" fmt, spaces_string__.c_str(), ##__VA_ARGS__);}


namespace Metris{
class DepthCounter{
  static int depth;
  static char spaces[65];
public:
  DepthCounter() = delete;

  DepthCounter(bool inc_) : inc(inc_){
    if(depth == -1) incDepth(); // initializes
    if(inc) incDepth();
  }
  ~DepthCounter(){
    if(inc) decDepth();
  }
  void decDepth(){ // can be called manually, for some rare cases like adaptGeoLines
    depth--;
    if(depth < 0) depth = 0;
    spaces[2*depth] = '\0';
    if(2*(depth + 1) < 64) spaces[2*(depth + 1)] = ' ';
  }
  void incDepth(){ // can be called manually, for some rare cases like adaptGeoLines
    depth++;
    METRIS_ASSERT(2*depth < 64);
    spaces[2*depth] = '\0';
    if(depth >= 1) spaces[2*(depth - 1)] = ' ';
  }
  //std::string getSpaces(){
  //  return " ";
  //}
  char* getSpaces(){
    return spaces;
  }
  int getDepth(){
    return depth;
  }

  void debug(){
    printf("-- debug dc depth = %d \n",depth);
    for(int ii = 0; ii < 65; ii++) printf(" char %d : -%c-\n",ii,spaces[ii]);
  }

  bool inc;
};


}//namespace
#endif