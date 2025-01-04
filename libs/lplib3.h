

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V3.73                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       Handles threads, scheduling, pipelines & dependencies */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 25 2008                                           */
/*   Last modification: nov 23 2021                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Types definitions                                                          */
/*----------------------------------------------------------------------------*/

#include <stdint.h>

#ifdef INT64
#define itg int64_t
#else
#define itg int32_t
#endif


/*----------------------------------------------------------------------------*/
/* User available procedures' prototypes                                      */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

int      AddDependency           (int64_t, itg, itg);
void     AddDependencyFast       (int64_t, int, itg *, int, itg *);
int      BeginDependency         (int64_t, int, int);
int      EndDependency           (int64_t, float [2]);
void     FreeType                (int64_t, int);
void     GetDependencyStats      (int64_t, int, int, float [2]);
void     GetLplibInformation     (int64_t, int *, int *);
int      GetNumberOfCores        ();
double   GetWallClock            ();
int      HilbertRenumbering      (int64_t, itg, double [6], double *, uint64_t *);
//int      HilbertRenumbering2D    (int64_t, itg, double [4], double *, uint64_t *);
int64_t  InitParallel            (int);
int64_t  InitParallelAttr        (int, size_t);
float    LaunchParallel          (int64_t, int, int, void *, void *);
float    LaunchParallelMultiArg  (int64_t, int, int, void *, int, ...);
int      LaunchPipeline          (int64_t, void *, void *, int, int *);
int      LaunchPipelineMultiArg  (int64_t, int, int *, void *prc, int, ...);
int      NewType                 (int64_t, itg);
int      ParallelMemClear        (int64_t, void *, size_t);
void     ParallelQsort           (int64_t, void *, size_t, size_t, int (*)(const void *, const void *));
int      ResizeType              (int64_t, int, itg);
void     StopParallel            (int64_t);
int      UpdateDependency        (int64_t, int, int, itg, itg);
void     UpdateDependencyFast    (int64_t, int, int, itg *, int, int, itg *);
void     WaitPipeline            (int64_t);
int      GetBlkIdx               (int64_t, int, int);
int      ChkBlkDep               (int64_t, int, int, int);
int      SetExtendedAttributes   (int64_t , ...);


#ifdef __cplusplus
} // end extern "C"
#endif


/*----------------------------------------------------------------------------*/
/* Public defines                                                             */
/*----------------------------------------------------------------------------*/

#define MaxPth 256

enum ArgAtr {
   SetInterleavingFactor = 1,
   SetInterleavingSize,
   DisableInterleaving,
   EnableBlockSorting,
   DisableBlockSorting,
   StaticScheduling
};

