

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                               LPlib V3.73                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Description:       Handles threads, scheduling                           */
/*                      & dependencies                                        */
/*   Author:            Loic MARECHAL                                         */
/*   Creation date:     feb 25 2008                                           */
/*   Last modification: nov 25 2021                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "lplib3.h"

#ifdef _WIN32
#include <windows.h>
#include <sys/timeb.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#if defined(_WIN32) && !defined(NATIVE_WINDOWS)
#include "winpthreads.h"
#else
#include <pthread.h>
#endif


/*----------------------------------------------------------------------------*/
/* Defines                                                                    */
/*----------------------------------------------------------------------------*/

#define MaxLibPar 10
#define MaxTyp    100
#define NmbSmlBlk 128
#define NmbDepBlk 512
#define MaxTotPip 65536
#define MaxPipDep 100
#define MaxHsh    10
#define HshBit    16
#define MaxF77Arg 20
#define WrkPerGrp 8

enum ParCmd {RunBigWrk, RunSmlWrk, RunDetWrk, ClrMem, EndPth};


/*----------------------------------------------------------------------------*/
/* Includes                                                                   */
/*----------------------------------------------------------------------------*/

#include <math.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>


/*----------------------------------------------------------------------------*/
/* Structures' prototypes                                                     */
/*----------------------------------------------------------------------------*/

typedef struct WrkSct
{
   itg               BegIdx, EndIdx, ItlTab[ MaxPth ][2];
   int               NmbDep, *DepWrdTab, GrpIdx;
   struct WrkSct     *pre, *nex;
}WrkSct;

typedef struct GrpSct
{
   int               idx, NmbSmlWrk[ MaxPth ];
   WrkSct            *SmlWrkTab[ MaxPth ][ WrkPerGrp ];
   struct GrpSct     *nex;
}GrpSct;

typedef struct
{
   itg               NmbLin, MaxNmbLin;
   int               NmbSmlWrk, SmlWrkSiz, DepWrkSiz, NmbGrp;
   int               NmbDepWrd, *DepWrdMat, *RunDepTab;
   WrkSct            *SmlWrkTab, *BigWrkTab;
   GrpSct            *NexGrp;
}TypSct;

typedef struct
{
   int               idx, NmbDetWrk;
   char              *ClrAdr;
   size_t            StkSiz;
   void *            *UsrStk;
   WrkSct            *wrk, **DetWrkTab;
   pthread_mutex_t   mtx;
   pthread_cond_t    cnd;
   pthread_t         pth;
   pthread_attr_t    atr;
   struct ParSct     *par;
}PthSct;

typedef struct PipSct
{
   int               idx, NmbF77Arg, NmbDep, DepTab[ MaxPipDep ];
   void              *prc, *arg, *F77ArgTab[ MaxF77Arg ];
   size_t            StkSiz;
   void              *UsrStk;
   pthread_attr_t    atr;
   pthread_t         pth;
   struct ParSct     *par;
}PipSct;

typedef struct ParSct
{
   int               NmbCpu, WrkCpt, NmbPip, PenPip, RunPip, NmbTyp, DynSch;
   int               req, cmd, *PipWrd, SizMul, NmbF77Arg, NmbVarArg;
   int               WrkSizSrt, NmbItlBlk, ItlBlkSiz, BufMax, BufCpt;
   size_t            StkSiz, ClrLinSiz;
   void              *F77ArgTab[ MaxF77Arg ];
   float             sta[2];
   void              (*prc)(itg, itg, int, void *), *arg;
   pthread_cond_t    ParCnd, PipCnd;
   pthread_mutex_t   ParMtx, PipMtx;
   pthread_t         PipPth;
   PthSct            *PthTab;
   TypSct            *TypTab, *CurTyp, *DepTyp, *typ1, *typ2;
   WrkSct            *NexWrk, *BufWrk[ MaxPth / 4 ];
}ParSct;

typedef struct
{
   uint64_t          *idx;
   double            box[6], *crd, *crd2;
}ArgSct;

typedef struct
{
   void              *base;
   size_t            nel, width;
   int               (*compar)(const void *, const void *);
}PipArgSct;


/*----------------------------------------------------------------------------*/
/* Private procedures' prototypes                                             */
/*----------------------------------------------------------------------------*/

static int     SetBit      (int *, int);
static int     GetBit      (int *, int);
static int     AndWrd      (int, int *, int *);
static void    AddWrd      (int, int *, int *);
static void    SubWrd      (int, int *, int *);
static void    ClrWrd      (int, int *);
static void    CpyWrd      (int, int *, int *);
int            CmpWrk      (const void *, const void *);
static void   *PipHdl      (void *);
static void   *PthHdl      (void *);
static WrkSct *NexWrk      (ParSct *, int);
void           PipSrt      (PipArgSct *);
static void    CalF77Prc   (itg, itg, int, ParSct *);
static void    CalF77Pip   (PipSct *, void *);
static void    CalVarArgPrc(itg, itg, int, ParSct *);
static int64_t IniPar      (int, size_t );
static void    SetItlBlk   (ParSct *, TypSct *);
static void    SetGrp      (ParSct *, TypSct *);


/*----------------------------------------------------------------------------*/
/* Init structures, scheduler and launch threads                              */
/*----------------------------------------------------------------------------*/

int64_t InitParallel(int NmbCpu)
{
   return(IniPar(NmbCpu, 0));
}

int64_t InitParallelAttr(int NmbCpu, size_t StkSiz)
{
#ifdef PTHREAD_STACK_MIN
   if(StkSiz < PTHREAD_STACK_MIN)
      StkSiz = PTHREAD_STACK_MIN;
#endif
   return(IniPar(NmbCpu, StkSiz));
}

static int64_t IniPar(int NmbCpu, size_t StkSiz)
{
   int i;
   int64_t ParIdx;
   ParSct *par;
   PthSct *pth;

   // Check the number of requested cpus
   if(NmbCpu < 1)
   {
      NmbCpu = GetNumberOfCores();

      if(NmbCpu < 1)
         return(0);
   }

   if(NmbCpu > MaxPth)
      NmbCpu = MaxPth;

   // Allocate and build main parallel structure
   if(!(par = calloc(1, sizeof(ParSct))))
      return(0);

   if(!(par->PthTab = calloc(NmbCpu, sizeof(PthSct))))
      return(0);

   if(!(par->TypTab = calloc((MaxTyp + 1), sizeof(TypSct))))
      return(0);

   if(!(par->PipWrd = calloc(MaxTotPip/32, sizeof(int))))
      return(0);

   par->NmbCpu = NmbCpu;
   par->WrkCpt = par->NmbPip = par->PenPip = par->RunPip = 0;
   par->SizMul = 2;
   par->StkSiz = StkSiz;
   par->NmbItlBlk = 1;
   par->WrkSizSrt = 1;
   par->DynSch = 1;

   // Set the size of WP buffer
   if(NmbCpu >= 4)
      par->BufMax = NmbCpu / 4;
   else
      par->BufMax = 1;

   pthread_mutex_init(&par->ParMtx, NULL);
   pthread_mutex_init(&par->PipMtx, NULL);
   pthread_cond_init(&par->ParCnd, NULL);
   pthread_cond_init(&par->PipCnd, NULL);

   // Launch pthreads
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pth->idx = i;
      pth->par = par;
      pthread_mutex_init(&pth->mtx, NULL);
      pthread_cond_init(&pth->cnd, NULL);

      if(StkSiz)
      {
         pthread_attr_init(&pth->atr);
         pth->StkSiz = StkSiz;
         pth->UsrStk = malloc(pth->StkSiz);
         pthread_attr_setstackaddr(&pth->atr, pth->UsrStk);
         pthread_attr_setstacksize(&pth->atr, pth->StkSiz);
         pthread_create(&pth->pth, &pth->atr, PthHdl, (void *)pth);
      }
      else
      {
         pth->StkSiz = 0;
         pth->UsrStk = NULL;
         pthread_create(&pth->pth, NULL, PthHdl, (void *)pth);
      }
   }

   // Wait for all threads to be up and wainting
   pthread_mutex_lock(&par->ParMtx);

   while(par->WrkCpt < par->NmbCpu)
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

   pthread_mutex_unlock(&par->ParMtx);

   ParIdx = (int64_t)par;

   return(ParIdx);
}


/*----------------------------------------------------------------------------*/
/* Stop all threads and free memories                                         */
/*----------------------------------------------------------------------------*/

void StopParallel(int64_t ParIdx)
{
   int i;
   PthSct *pth;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Send stop to all threads
   pthread_mutex_lock(&par->ParMtx);
   par->cmd = EndPth;
   pthread_mutex_unlock(&par->ParMtx);

   // Wait for all threads to complete
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pthread_mutex_lock(&pth->mtx);
      pthread_cond_signal(&pth->cnd);
      pthread_mutex_unlock(&pth->mtx);
      pthread_join(pth->pth, NULL);

      if(pth->UsrStk)
         free(pth->UsrStk);
   }

   pthread_mutex_destroy(&par->ParMtx);
   pthread_cond_destroy(&par->ParCnd);

   WaitPipeline(ParIdx);

   pthread_mutex_destroy(&par->PipMtx);
   pthread_cond_destroy(&par->PipCnd);

   // Free memories
   for(i=1;i<=MaxTyp;i++)
      if(par->TypTab[i].NmbLin)
         FreeType(ParIdx, i);

   free(par->PthTab);
   free(par->TypTab);
   free(par->PipWrd);
   free(par);
}


/*----------------------------------------------------------------------------*/
/* Try to detect the number of system's threads                               */
/*----------------------------------------------------------------------------*/

int GetNumberOfCores()
{
#ifdef _WIN32
   SYSTEM_INFO info;
   GetSystemInfo(&info);
   return(info.dwNumberOfProcessors);
#else
#ifdef _SC_NPROCESSORS_ONLN
   return((itg)sysconf(_SC_NPROCESSORS_ONLN));
#else
   return(1);
#endif
#endif
}


/*----------------------------------------------------------------------------*/
/* Set extra attributes to fine tune the blocks processing order              */
/*----------------------------------------------------------------------------*/

int SetExtendedAttributes(int64_t ParIdx, ...)
{
   int NmbArg = 0, ArgCod, ArgVal;
   ParSct *par = (ParSct *)ParIdx;
   va_list ArgLst;

   // Attributes cannot be modified the a prallel loop is running
   if(par->typ1)
      return(0);

   // Read the argument list that must be terminated by a zero
   va_start(ArgLst, ParIdx);

   do
   {
      ArgCod = va_arg(ArgLst, int);

      switch(ArgCod)
      {
         // Set the number of interleaved blocks in independant loops
         case SetInterleavingFactor :
         {
            ArgVal = va_arg(ArgLst, int);

            if(ArgVal > 0)
            {
               par->NmbItlBlk = ArgVal;
               par->ItlBlkSiz = 0;
               NmbArg++;
            }
         }break;

         // Set the interleaved blocks size in independant loops
         case SetInterleavingSize :
         {
            ArgVal = va_arg(ArgLst, int);

            if(ArgVal > 0)
            {
               par->NmbItlBlk = 0;
               par->ItlBlkSiz = ArgVal;
               NmbArg++;
            }
         }break;

         // Desable blocks interleaving (default)
         case DisableInterleaving :
         {
            par->NmbItlBlk = 1;
            par->ItlBlkSiz = 0;
            NmbArg++;
         }break;

         // Sort depedancy-loop WP through their number of dependencies (default)
         case EnableBlockSorting :
         {
            par->WrkSizSrt = 1;
            NmbArg++;
         }break;

         // Disable WP sorting: it lowers concurrency but enhances cache reuse
         case DisableBlockSorting :
         {
            par->WrkSizSrt = 0;
            NmbArg++;
         }break;

         // Static scheduling makes the library deterministic
         case StaticScheduling :
         {
            // WP sorting is useless in this mode so it is disabled
            par->WrkSizSrt = par->DynSch = 0;
            NmbArg++;
         }break;
      }
   }while(ArgCod);

   va_end(ArgLst);

   return(NmbArg);
}


/*----------------------------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2                      */
/*----------------------------------------------------------------------------*/

float LaunchParallel(int64_t ParIdx, int TypIdx1, int TypIdx2,
                     void *prc, void *PtrArg )
{
   int      i;
   float    acc;
   PthSct   *pth;
   ParSct   *par = (ParSct *)ParIdx;
   TypSct   *typ1, *typ2 = NULL;
   GrpSct   *grp;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1.);

   // Check bounds
   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 0)
   ||  (TypIdx2 > MaxTyp) || (TypIdx1 == TypIdx2) )
   {
      return(-1.);
   }

   typ1 =  &par->TypTab[ TypIdx1 ];

   // Launch small WP with static scheduling
   if(TypIdx2 && !par->DynSch)
   {
      grp = typ1->NexGrp;
      acc = 0.;

      do
      {
         // Lock acces to global parameters
         pthread_mutex_lock(&par->ParMtx);

         par->cmd = RunDetWrk;
         par->prc = (void (*)(itg, itg, int, void *))prc;
         par->arg = PtrArg;
         par->typ1 = typ1;
         par->typ2 = NULL;
         par->WrkCpt = 0;

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];
            pth->NmbDetWrk = grp->NmbSmlWrk[i];
            pth->DetWrkTab = grp->SmlWrkTab[i];
            acc += (float)grp->NmbSmlWrk[i];
         }

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];
            pthread_mutex_lock(&pth->mtx);
            pthread_cond_signal(&pth->cnd);
            pthread_mutex_unlock(&pth->mtx);
         }

         pthread_cond_wait(&par->ParCnd, &par->ParMtx);

         pthread_mutex_unlock(&par->ParMtx);
         grp = grp->nex;
      }while(grp);

      acc /= (float)(NmbSmlBlk * typ1->NmbGrp) / (float)WrkPerGrp;
   }
   else if(TypIdx2 && par->DynSch)
   {
      // Launch small WP with dynamic scheduling

      // Lock acces to global parameters
      pthread_mutex_lock(&par->ParMtx);

      par->cmd = RunSmlWrk;
      par->prc = (void (*)(itg, itg, int, void *))prc;
      par->arg = PtrArg;
      par->typ1 = typ1;
      par->typ2 = typ2 = &par->TypTab[ TypIdx2 ];
      par->NexWrk = typ1->SmlWrkTab;
      par->BufCpt = 0;
      par->WrkCpt = 0;
      par->sta[0] = par->sta[1] = 0.;
      par->req = 0;

      // Clear running wp
      for(i=0;i<par->NmbCpu;i++)
         par->PthTab[i].wrk = NULL;

      ClrWrd(typ1->NmbDepWrd, typ1->RunDepTab);

      // Build a linked list of wp
      for(i=0;i<par->typ1->NmbSmlWrk;i++)
      {
         typ1->SmlWrkTab[i].pre = &typ1->SmlWrkTab[ i-1 ];
         typ1->SmlWrkTab[i].nex = &typ1->SmlWrkTab[ i+1 ];
      }

      typ1->SmlWrkTab[0].pre = typ1->SmlWrkTab[ typ1->NmbSmlWrk - 1 ].nex = NULL;

      // Main loop: wake up threads and wait for completion or blocked threads
      do
      {
         // Search for some idle threads
         par->req = 0;

         for(i=0;i<par->NmbCpu;i++)
         {
            pth = &par->PthTab[i];

            if(pth->wrk)
               continue;

            if(!(pth->wrk = NexWrk(par, i)))
            {
               par->req = 1;
               break;
            }

            // Wake up the thread and provide it with a WP list
            pthread_mutex_lock(&pth->mtx);
            pthread_cond_signal(&pth->cnd);
            pthread_mutex_unlock(&pth->mtx);
         }

         // If every WP are done : exit the parallel loop
         if(par->WrkCpt == typ1->NmbSmlWrk)
            break;

         // Otherwise, wait for a blocked thread
         pthread_cond_wait(&par->ParCnd, &par->ParMtx);
      }while(1);

      pthread_mutex_unlock(&par->ParMtx);

      // Compute the average concurrency factor
      acc = par->sta[0] ? (par->sta[1] / par->sta[0]) : 0;
   }
   else
   {
      // Launch big WP with static scheduling

      // Lock acces to global parameters
      pthread_mutex_lock(&par->ParMtx);

      par->cmd = RunBigWrk;
      par->prc = (void (*)(itg, itg, int, void *))prc;
      par->arg = PtrArg;
      par->typ1 = typ1;
      par->typ2 = NULL;
      par->WrkCpt = 0;

      for(i=0;i<par->NmbCpu;i++)
      {
         pth = &par->PthTab[i];
         pth->wrk = &typ1->BigWrkTab[i];
      }

      // Update block interleaving according to the current attributes
      SetItlBlk(par, typ1);

      for(i=0;i<par->NmbCpu;i++)
      {
         pth = &par->PthTab[i];
         pthread_mutex_lock(&pth->mtx);
         pthread_cond_signal(&pth->cnd);
         pthread_mutex_unlock(&pth->mtx);
      }

      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

      pthread_mutex_unlock(&par->ParMtx);

      // Arbitrary set the average concurrency factor
      acc = (float)par->NmbCpu;
   }

   // Clear the main datatyp loop to indicate that no LaunchParallel is running
   par->typ1 = 0;

   // Return the concurrency factor
   return(acc);
}


/*----------------------------------------------------------------------------*/
/* Launch a parallel procudure with variable arguments.                       */
/* Arguments are passed as pointer to void.                                   */
/*----------------------------------------------------------------------------*/

float LaunchParallelMultiArg( int64_t ParIdx, int TypIdx1, int TypIdx2,
                              void *prc, int NmbArg, ... )
{
   int i;
   float acc;
   va_list ArgLst;
   ParSct *par = (ParSct *)ParIdx;

   if(NmbArg > 20)
      return(-1.);

   par->NmbVarArg = NmbArg;
   va_start(ArgLst, NmbArg);

   for(i=0;i<NmbArg;i++)
      par->F77ArgTab[i] = va_arg(ArgLst, void *);

   va_end(ArgLst);

   acc = LaunchParallel(ParIdx, TypIdx1, TypIdx2, prc, NULL);
   par->NmbVarArg = 0;

   return(acc);
}


/*----------------------------------------------------------------------------*/
/* Pthread handler, waits for job, does it, then signal end                   */
/*----------------------------------------------------------------------------*/

static void *PthHdl(void *ptr)
{
   itg i, beg, end;
   PthSct *pth = (PthSct *)ptr;
   ParSct *par = pth->par;

   // Tell the scheduler if all threads are ready
   pthread_mutex_lock(&par->ParMtx);
   par->WrkCpt++;
   pthread_cond_signal(&par->ParCnd);
   pthread_mutex_lock(&pth->mtx);
   pthread_mutex_unlock(&par->ParMtx);

   // Enter main loop until StopParallel is send
   do
   {
      // Wait for a wake-up signal from the main loop
      pthread_cond_wait(&pth->cnd, &pth->mtx);

      switch(par->cmd)
      {
         // Call user's procedure with big WP
         case RunBigWrk :
         {
            // Loop over the interleaved blocks
            for(i=0;i<par->NmbItlBlk;i++)
            {
               beg = pth->wrk->ItlTab[i][0];
               end = pth->wrk->ItlTab[i][1];

               if(!beg || !end || (end < beg))
                  continue;

               // Launch a single big wp and signal completion to the scheduler
               if(par->NmbF77Arg)
                  CalF77Prc(beg, end, pth->idx, par);
               else if(par->NmbVarArg)
                  CalVarArgPrc(beg, end, pth->idx, par);
               else
                  par->prc(beg, end, pth->idx, par->arg);
            }

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;

            if(par->WrkCpt >= par->NmbCpu)
               pthread_cond_signal(&par->ParCnd);

            pthread_mutex_unlock(&par->ParMtx);
         }break;

         // Call user's procedure with small WP using dynamic scheduling
         case RunSmlWrk :
         {
            do
            {
               // Run the WP
               if(par->NmbF77Arg)
                  CalF77Prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par);
               else if(par->NmbVarArg)
                  CalVarArgPrc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par);
               else
                  par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

               // Locked acces to global parameters: 
               // update WP count, tag WP done and signal the main loop
               pthread_mutex_lock(&par->ParMtx);

               par->WrkCpt++;

               if(!(pth->wrk = NexWrk(par, pth->idx)))
               {
                  par->req = 1;
                  pthread_cond_signal(&par->ParCnd);
                  pthread_mutex_unlock(&par->ParMtx);
                  break;
               }

               if(par->req)
                  pthread_cond_signal(&par->ParCnd);

               pthread_mutex_unlock(&par->ParMtx);
            }while(1);
         }break;

         // Call user's procedure with small WP using static scheduling
         case RunDetWrk :
         {
            // Loop over the groups' WP
            for(i=0;i<pth->NmbDetWrk;i++)
            {
               beg = pth->DetWrkTab[i]->BegIdx;
               end = pth->DetWrkTab[i]->EndIdx;

               if(par->NmbF77Arg)
                  CalF77Prc(beg, end, pth->idx, par);
               else if(par->NmbVarArg)
                  CalVarArgPrc(beg, end, pth->idx, par);
               else
                  par->prc(beg, end, pth->idx, par->arg);
            }

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;

            if(par->WrkCpt >= par->NmbCpu)
               pthread_cond_signal(&par->ParCnd);

            pthread_mutex_unlock(&par->ParMtx);
         }break;

         case ClrMem :
         {
            // Clear memory and signal completion to the scheduler
            memset(pth->ClrAdr, 0, par->ClrLinSiz);

            pthread_mutex_lock(&par->ParMtx);
            par->WrkCpt++;
            pthread_cond_signal(&par->ParCnd);
            pthread_mutex_unlock(&par->ParMtx);
         }break;

         case EndPth :
         {
            // Destroy the thread mutex and condition and call for join
            pthread_mutex_unlock(&pth->mtx);
            pthread_mutex_destroy(&pth->mtx);
            pthread_cond_destroy(&pth->cnd);
            return(NULL);
         }break;
      }
   }while(1);

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Get the next WP to be computed                                             */
/*----------------------------------------------------------------------------*/

static WrkSct *NexWrk(ParSct *par, int PthIdx)
{
   int i;
   PthSct *pth = &par->PthTab[ PthIdx ];
   WrkSct *wrk;

   // Update stats
   par->sta[0]++;

   for(i=0;i<par->NmbCpu;i++)
      if(par->PthTab[i].wrk)
         par->sta[1]++;

   // Remove previous work's tags
   if(pth->wrk)
      SubWrd(par->typ1->NmbDepWrd, pth->wrk->DepWrdTab, par->typ1->RunDepTab);

   // If the wp's buffer is empty search for some new compatible wp to fill in
   if(!par->BufCpt)
   {
      wrk = par->NexWrk;

      while(wrk)
      {
         // Check for dependencies
         if(!AndWrd(par->typ1->NmbDepWrd, wrk->DepWrdTab, par->typ1->RunDepTab))
         {
            par->BufWrk[ par->BufCpt++ ] = wrk;

            // Unlink wp
            if(wrk->pre)
               wrk->pre->nex = wrk->nex;
            else
               par->NexWrk = wrk->nex;

            if(wrk->nex)
               wrk->nex->pre = wrk->pre;

            // Add new work's tags
            AddWrd(par->typ1->NmbDepWrd, wrk->DepWrdTab, par->typ1->RunDepTab);

            if(par->BufCpt == par->BufMax)
               break;
         }

         wrk = wrk->nex;
      }
   }

   // Return the next available wp in buffer and unlink it from the todo list
   return(par->BufCpt ? par->BufWrk[ --par->BufCpt ] : NULL);
}


/*----------------------------------------------------------------------------*/
/* Allocate a new kind of elements and set work-packages                      */
/*----------------------------------------------------------------------------*/

int NewType(int64_t ParIdx, itg NmbLin)
{
   int      TypIdx = 0;
   itg      i, idx;
   TypSct   *typ;
   ParSct   *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   if(NmbLin <= 0)
      return(0);

   // Search for a free type structure
   for(i=1;i<=MaxTyp;i++)
      if(!par->TypTab[i].NmbLin)
      {
         TypIdx = i;
         break;
      }

   if(!TypIdx)
      return(0);

   typ = &par->TypTab[ TypIdx ];
   typ->NmbLin = NmbLin;
   typ->MaxNmbLin = NmbLin * par->SizMul;

   // Compute the size of small work-packages
   if(NmbLin >= NmbSmlBlk * par->NmbCpu)
   {
      typ->SmlWrkSiz = NmbLin / (NmbSmlBlk * par->NmbCpu);
      typ->NmbSmlWrk = NmbLin / typ->SmlWrkSiz;

      if(NmbLin != typ->NmbSmlWrk * typ->SmlWrkSiz)
         typ->NmbSmlWrk++;
   }
   else
   {
      typ->SmlWrkSiz = NmbLin;
      typ->NmbSmlWrk = 1;
   }

   if(!(typ->SmlWrkTab = calloc(typ->NmbSmlWrk * par->SizMul , sizeof(WrkSct))))
      return(0);

   // Set small work-packages
   idx = 0;

   for(i=0;i<typ->NmbSmlWrk;i++)
   {
      typ->SmlWrkTab[i].BegIdx = idx + 1;
      typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
      idx += typ->SmlWrkSiz;
   }

   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

   // Compute the size of big work-packages
   if(!(typ->BigWrkTab = calloc(par->NmbCpu * par->SizMul , sizeof(WrkSct))))
      return(0);

   return(TypIdx);
}


/*----------------------------------------------------------------------------*/
/* Resize a data type up to a two fold increase                               */
/*----------------------------------------------------------------------------*/

int ResizeType(int64_t ParIdx, int TypIdx, itg NmbLin)
{
   itg      i, idx;
   TypSct   *typ;
   ParSct   *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds and free mem
   if( (TypIdx < 1) || (TypIdx > MaxTyp) )
      return(0);

   typ = &par->TypTab[ TypIdx ];

   if( (NmbLin < typ->NmbLin) || (NmbLin > typ->MaxNmbLin) )
      return(0);

   // Set small work-packages
   i = typ->NmbSmlWrk;
   idx = typ->NmbLin;
   typ->NmbLin = NmbLin;

   while(idx < NmbLin)
   {
      typ->SmlWrkTab[i].BegIdx = idx + 1;
      typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
      i++;
      idx += typ->SmlWrkSiz;
      typ->NmbSmlWrk++;
   }

   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

   return(TypIdx);
}


/*----------------------------------------------------------------------------*/
/* Set big work blocks interleaving                                           */
/*----------------------------------------------------------------------------*/

static void SetItlBlk(ParSct *par, TypSct *typ)
{
   itg      i, j, BegIdx, EndIdx, CpuIdx = 0, PagIdx[ MaxPth ] = {0};
   double   ItlSiz, ItlIdx = 0.;

   // Set big WP interleaved indices and block sizes if requested
   if(par->NmbItlBlk)
      ItlSiz = (double)typ->NmbLin / (double)(par->NmbItlBlk * par->NmbCpu);
   else if(par->ItlBlkSiz)
   {
      par->NmbItlBlk = typ->NmbLin / (par->ItlBlkSiz * par->NmbCpu);
      ItlSiz = (double)par->ItlBlkSiz;
   }

   // In case the block or lines numbers is too small, deactivate interleaving
   if(!par->NmbItlBlk || !ItlSiz)
   {
      par->NmbItlBlk = 1;
      ItlSiz = (double)typ->NmbLin / (double)(par->NmbCpu);
   }

   for(i=0;i<par->NmbCpu;i++)
   {
      for(j=0;j<par->NmbItlBlk;j++)
      {
         BegIdx = (itg)(ItlIdx + 1.);
         EndIdx = (itg)(ItlIdx + ItlSiz);
         ItlIdx += ItlSiz;

         if(BegIdx <= EndIdx)
         {
            typ->BigWrkTab[ CpuIdx ].ItlTab[ PagIdx[ CpuIdx ] ][0] = BegIdx;
            typ->BigWrkTab[ CpuIdx ].ItlTab[ PagIdx[ CpuIdx ] ][1] = EndIdx;
            PagIdx[ CpuIdx ]++;
            CpuIdx++;
            CpuIdx = CpuIdx % par->NmbCpu;
         }
      }
   }

   typ->BigWrkTab[ par->NmbCpu - 1 ].ItlTab[ par->NmbItlBlk - 1 ][1] = typ->NmbLin;
}


/*----------------------------------------------------------------------------*/
/* Add this kind of element to the free-list                                  */
/*----------------------------------------------------------------------------*/

void FreeType(int64_t ParIdx, int TypIdx)
{
   TypSct *typ;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Check bounds and free mem
   if( (TypIdx < 1) || (TypIdx > MaxTyp) )
      return;

   typ = &par->TypTab[ TypIdx ];

   if(typ->SmlWrkTab)
      free(typ->SmlWrkTab);

   if(typ->BigWrkTab)
      free(typ->BigWrkTab);

   if(typ->RunDepTab)
      free(typ->RunDepTab);

   if(typ->DepWrdMat)
      free(typ->DepWrdMat);

   memset(typ, 0, sizeof(TypSct));
}


/*----------------------------------------------------------------------------*/
/* Allocate a dependency matrix linking both types                            */
/*----------------------------------------------------------------------------*/

int BeginDependency(int64_t ParIdx, int TypIdx1, int TypIdx2)
{
   int i;
   TypSct *typ1, *typ2;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   par->CurTyp = typ1 = &par->TypTab[ TypIdx1 ];
   par->DepTyp = typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin || !typ2->NmbLin)
   {
      return(0);
   }

   // Compute dependency table's size
   if( (typ2->NmbLin >= NmbDepBlk * par->NmbCpu)
   &&  (typ2->NmbLin >= typ1->DepWrkSiz * 32) )
   {
      typ1->DepWrkSiz = typ2->NmbLin / (NmbDepBlk * par->NmbCpu);
      typ1->NmbDepWrd = typ2->NmbLin / (typ1->DepWrkSiz * 32);

      if(typ2->NmbLin != typ1->NmbDepWrd * typ1->DepWrkSiz * 32)
         typ1->NmbDepWrd++;
   }
   else
   {
      typ1->DepWrkSiz = typ2->NmbLin;
      typ1->NmbDepWrd = 1;
   }

   // Allocate a global dependency table
   if(!(typ1->DepWrdMat =
      calloc(typ1->NmbSmlWrk * typ1->NmbDepWrd * par->SizMul, sizeof(int))))
   {
      return(0);
   }

   // Then spread sub-tables among WP
   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      typ1->SmlWrkTab[i].NmbDep = 0;
      typ1->SmlWrkTab[i].DepWrdTab =
         &typ1->DepWrdMat[ i * typ1->NmbDepWrd * par->SizMul ];
   }

   // Allocate a running tags table
   if(!(typ1->RunDepTab = calloc(typ1->NmbDepWrd * par->SizMul, sizeof(int))))
      return(0);

   return(typ1->NmbDepWrd);
}


/*----------------------------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2                           */
/*----------------------------------------------------------------------------*/

int AddDependency(int64_t ParIdx, itg idx1, itg idx2)
{
   WrkSct *wrk;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if( !par || !par->CurTyp || !par->CurTyp->SmlWrkSiz
   ||  !par->CurTyp->DepWrkSiz || !par->DepTyp || (idx1 < 1)
   ||  (idx1 > par->CurTyp->NmbLin) || (idx2 < 1)
   ||  (idx2 > par->DepTyp->NmbLin) )
   {
      return(0);
   }

   // Set and count dependency bit
   wrk = &par->CurTyp->SmlWrkTab[ (idx1-1) / par->CurTyp->SmlWrkSiz ];

   if(!SetBit(wrk->DepWrdTab, (idx2-1) / par->CurTyp->DepWrkSiz ))
      wrk->NmbDep++;

   return(wrk->NmbDep);
}


/*----------------------------------------------------------------------------*/
/* Set all to all dependencies without any pre-checking                       */
/*----------------------------------------------------------------------------*/

void AddDependencyFast( int64_t ParIdx, int NmbTyp1, itg *TabIdx1,
                        int NmbTyp2, itg *TabIdx2 )
{
   int i, j;
   ParSct *par = (ParSct *)ParIdx;
   WrkSct *wrk;

   for(i=0;i<NmbTyp1;i++)
   {
      wrk = &par->CurTyp->SmlWrkTab[ (TabIdx1[i] - 1) / par->CurTyp->SmlWrkSiz ];

      for(j=0;j<NmbTyp2;j++)
         if( !SetBit(wrk->DepWrdTab, (TabIdx2[j] - 1) / par->CurTyp->DepWrkSiz ) )
            wrk->NmbDep++;
   }
}


/*----------------------------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2                           */
/*----------------------------------------------------------------------------*/

int UpdateDependency(int64_t ParIdx, int TypIdx1, int TypIdx2,
                     itg idx1, itg idx2 )
{
   WrkSct *wrk;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(0);

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin || !typ2->NmbLin
   ||  !typ1->SmlWrkSiz || !typ1->DepWrkSiz || (idx1 < 1)
   ||  (idx1 > typ1->MaxNmbLin) || (idx2 < 1) || (idx2 > typ2->MaxNmbLin) )
   {
      return(0);
   }

   // Set and count dependency bit
   wrk = &typ1->SmlWrkTab[ (idx1-1) / typ1->SmlWrkSiz ];

   if(!SetBit(wrk->DepWrdTab, (idx2-1) / typ1->DepWrkSiz ))
      wrk->NmbDep++;

   return(wrk->NmbDep);
}


/*----------------------------------------------------------------------------*/
/* Same as above, without any pre-checking of data                            */
/*----------------------------------------------------------------------------*/

void UpdateDependencyFast( int64_t ParIdx,  int TypIdx1, int NmbTyp1,
                           itg *TabIdx1, int TypIdx2, int NmbTyp2,
                           itg *TabIdx2 )
{
   int i, j;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1 = &par->TypTab[ TypIdx1 ];
   WrkSct *wrk;
   (void)(TypIdx2);

   for(i=0;i<NmbTyp1;i++)
   {
      wrk = &typ1->SmlWrkTab[ (TabIdx1[i] - 1) / typ1->SmlWrkSiz ];

      for(j=0;j<NmbTyp2;j++)
         if( !SetBit(wrk->DepWrdTab, (TabIdx2[j] - 1) / typ1->DepWrkSiz ) )
            wrk->NmbDep++;
   }
}


/*----------------------------------------------------------------------------*/
/* Sort wp depending on their number of dependencies                          */
/*----------------------------------------------------------------------------*/

int EndDependency(int64_t ParIdx, float DepSta[2])
{
   int      i, NmbDepBit, TotNmbDep = 0;
   ParSct   *par = (ParSct *)ParIdx;
   TypSct   *typ1, *typ2;

   // Get and check lib parallel instance
   if( !ParIdx || !DepSta )
      return(0);

   // Compute average number of collisions
   DepSta[1] = 0.;
   typ1 = par->CurTyp;
   typ2 = par->DepTyp;

   if(!typ1 || !typ2 || !typ1->DepWrkSiz)
      return(0);

   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      TotNmbDep += typ1->SmlWrkTab[i].NmbDep;

      if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
         DepSta[1] = (float)typ1->SmlWrkTab[i].NmbDep;
   }

   if(!TotNmbDep)
      return(0);

   DepSta[0] = (float)TotNmbDep;

   // Compute stats
   if(typ2->NmbLin >= typ1->DepWrkSiz)
   {
      NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

      if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
         NmbDepBit++;
   }
   else
      NmbDepBit = 1;

   if(!NmbDepBit)
      return(0);

   DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
   DepSta[1] = 100 * DepSta[1] / NmbDepBit;

   // Sort WP from highest collision number to the lowest
   if(par->WrkSizSrt && par->DynSch)
      qsort(typ1->SmlWrkTab, typ1->NmbSmlWrk, sizeof(WrkSct), CmpWrk);

   if(!par->DynSch)
      SetGrp(par, typ1);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Sort wp depending on their number of dependencies                          */
/*----------------------------------------------------------------------------*/

void GetDependencyStats(int64_t ParIdx, int TypIdx1,
                        int TypIdx2, float DepSta[2])
{
   int i, NmbDepBit, TotNmbDep = 0;
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ1, *typ2;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   // Check bounds
   typ1 = &par->TypTab[ TypIdx1 ];
   typ2 = &par->TypTab[ TypIdx2 ];

   if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1)
   ||  (TypIdx2 > MaxTyp) || (typ1 == typ2) || !typ1->NmbLin
   ||  !typ2->NmbLin )
   {
      return;
   }

   // Compute average number of collisions
   DepSta[1] = 0.;

   for(i=0;i<typ1->NmbSmlWrk;i++)
   {
      TotNmbDep += typ1->SmlWrkTab[i].NmbDep;

      if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
         DepSta[1] = (float)typ1->SmlWrkTab[i].NmbDep;
   }

   DepSta[0] = (float)TotNmbDep;

   // Compute stats
   if(typ2->NmbLin >= typ1->DepWrkSiz)
   {
      NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

      if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
         NmbDepBit++;
   }
   else
      NmbDepBit = 1;

   DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
   DepSta[1] = 100 * DepSta[1] / NmbDepBit;
}


/*----------------------------------------------------------------------------*/
/* Return the block ID containing the given element ID                        */
/*----------------------------------------------------------------------------*/

int GetBlkIdx(int64_t ParIdx, int typ, int idx)
{
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1);

   return( (idx - 1) / par->TypTab[ typ ].SmlWrkSiz );
}


/*----------------------------------------------------------------------------*/
/* Check whether two blocks are dependant or not from each others             */
/*----------------------------------------------------------------------------*/

int ChkBlkDep(int64_t ParIdx, int TypIdx, int blk1, int blk2)
{
   ParSct *par = (ParSct *)ParIdx;
   TypSct *typ;

   // Get and check lib parallel instance
   if(!ParIdx)
      return(-1);

   typ = &par->TypTab[ TypIdx ];

   return(AndWrd( typ->NmbDepWrd,
                  typ->SmlWrkTab[ blk1 ].DepWrdTab,
                  typ->SmlWrkTab[ blk2 ].DepWrdTab ));
}


/*----------------------------------------------------------------------------*/
/* Test and set a bit in a multibyte word                                     */
/*----------------------------------------------------------------------------*/

static int SetBit(int *tab, int idx)
{
   int res = ( tab[ idx >> 5 ] & (1UL << (idx & 31)) );
   tab[ idx >> 5 ] |= 1UL << (idx & 31);
   return(res);
}


/*----------------------------------------------------------------------------*/
/* Test a bit in a multibyte word                                             */
/*----------------------------------------------------------------------------*/

static int GetBit(int *tab, int idx)
{
   return( tab[ idx >> 5 ] & (1UL << (idx & 31)) );
}


/*----------------------------------------------------------------------------*/
/* Check wether two WP share common resources -> locked                       */
/*----------------------------------------------------------------------------*/

static int AndWrd(int NmbWrd, int *wrd1, int *wrd2)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      if(wrd1[i] & wrd2[i])
         return(1);

   return(0);
}


/*----------------------------------------------------------------------------*/
/* Logical OR between two multibyte words                                     */
/*----------------------------------------------------------------------------*/

static void AddWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] |= SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Exclusive OR between two multibyte words                                   */
/*----------------------------------------------------------------------------*/

static void SubWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] &= ~SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Clear a multibyte word                                                     */
/*----------------------------------------------------------------------------*/

static void ClrWrd(int NmbWrd, int *wrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      wrd[i] = 0;
}


/*----------------------------------------------------------------------------*/
/* Copy a multibyte word into another one                                     */
/*----------------------------------------------------------------------------*/

static void CpyWrd(int NmbWrd, int *SrcWrd, int *DstWrd)
{
   int i;

   for(i=0;i<NmbWrd;i++)
      DstWrd[i] = SrcWrd[i];
}


/*----------------------------------------------------------------------------*/
/* Compare two workpackages number of bits                                    */
/*----------------------------------------------------------------------------*/

int CmpWrk(const void *ptr1, const void *ptr2)
{
   WrkSct *w1, *w2;

   w1 = (WrkSct *)ptr1;
   w2 = (WrkSct *)ptr2;

   if(w1->NmbDep > w2->NmbDep)
      return(-1);
   else if(w1->NmbDep < w2->NmbDep)
      return(1);
   else
      return(0);
}


/*----------------------------------------------------------------------------*/
/* Generate static scheduling groups of small WP for each thread              */
/*----------------------------------------------------------------------------*/

static void SetGrp(ParSct *par, TypSct *typ)
{
   int      i, NmbSmlWrk, *GrpWrd, *AllWrd, *TstWrd;
   int      IncFlg, siz = typ->NmbDepWrd;
   GrpSct   *grp;
   WrkSct   *wrk, *NexWrk;

   // Initialize the static group list
   NmbSmlWrk = typ->NmbSmlWrk;
   typ->NmbGrp = 0;
   typ->NexGrp = NULL;

   // Allocate a dependency word to contain all threads
   GrpWrd = malloc(par->NmbCpu * siz * sizeof(int));
   assert(GrpWrd);

   // Allocate a dependency word to concatenate all thread words
   AllWrd = malloc(siz * sizeof(int));
   assert(AllWrd);

   // Allocate a working dependency word for testings
   TstWrd = malloc(siz * sizeof(int));
   assert(TstWrd);

   // Link all WP together to make a free list
   NexWrk = &typ->SmlWrkTab[0];

   for(i=0;i<typ->NmbSmlWrk;i++)
   {
      typ->SmlWrkTab[i].pre = &typ->SmlWrkTab[ i - 1 ];
      typ->SmlWrkTab[i].nex = &typ->SmlWrkTab[ i + 1 ];
   }

   typ->SmlWrkTab[0].pre = NULL;
   typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].nex = NULL;

   // Kepp on creating new groups as long as there are free WP
   do
   {
      // Allocate a new group and link it
      grp = calloc(1, sizeof(GrpSct));
      assert(grp);
      grp->nex = typ->NexGrp;
      typ->NexGrp = grp;
      grp->idx = ++typ->NmbGrp;

      // Reset all the dependencies words
      memset(GrpWrd, 0, par->NmbCpu * siz * sizeof(int));
      memset(AllWrd, 0, siz * sizeof(int));

      // Keep on adding WP to each groups' threads until they reach WrkPerGrp
      do
      {
         // Could we add some WP to this group during this iteration ?
         IncFlg = 0;

         // Loop over the threads and try adding them as many WP as possible
         for(i=0;i<par->NmbCpu;i++)
         {
            // Do not add WP is this threads has reached its target number
            if(!NexWrk || (grp->NmbSmlWrk[i] >= WrkPerGrp))
               continue;

            // Initialize the WP pointer with the first
            // available one from the linked list
            wrk = NexWrk;

            // Then loop over the linked list of available WP and add them to
            // the curent group-thread list as long as they do not interfere
            // with the combined dependency word of other threads
            do
            {
               // To avoid making a AND with all threads' words, we remode (XOR)
               // the tested thread' word from the combined one and test (AND)
               // it against this WP's word
               CpyWrd(siz, AllWrd, TstWrd);
               SubWrd(siz, &GrpWrd[ i * siz ], TstWrd);

               if(!AndWrd(siz, wrk->DepWrdTab, TstWrd))
               {
                  // If this WP is compatible, add its word to the thread
                  // dependency word and to the combined one
                  AddWrd(siz, wrk->DepWrdTab, AllWrd);
                  AddWrd(siz, wrk->DepWrdTab, &GrpWrd[ i * siz ]);

                  // Add this WP to the list, decrease the number of available
                  // WPs and set the flag to indicate that we found some work
                  grp->SmlWrkTab[i][ grp->NmbSmlWrk[i] ] = wrk;
                  grp->NmbSmlWrk[i]++;
                  NmbSmlWrk--;
                  IncFlg = 1;

                  // Unlink the WP from the available WP linked list
                  if(wrk->pre)
                     wrk->pre->nex = wrk->nex;
                  else
                     NexWrk = wrk->nex;

                  if(wrk->nex)
                     wrk->nex->pre = wrk->pre;
               }

               // Stop adding WP to this group-thread if its table is full
               if(grp->NmbSmlWrk[i] >= WrkPerGrp)
                  break;

               wrk = wrk->nex;
            }while(wrk && NmbSmlWrk);
         }
      }while(IncFlg && NmbSmlWrk);
   }while(NmbSmlWrk);

   free(GrpWrd);
}


/*----------------------------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2                      */
/*----------------------------------------------------------------------------*/

int ParallelMemClear(int64_t ParIdx, void *PtrArg, size_t siz)
{
   char *tab = (char *)PtrArg;
   int i;
   PthSct *pth;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance, adresse and size
   if(!ParIdx || !tab || (siz < (size_t)par->NmbCpu) )
      return(0);

   // Lock acces to global parameters
   pthread_mutex_lock(&par->ParMtx);

   par->cmd = ClrMem;
   par->ClrLinSiz = siz / par->NmbCpu;
   par->WrkCpt = 0;

   // Spread the buffer among each thread and wake then up
   for(i=0;i<par->NmbCpu;i++)
   {
      pth = &par->PthTab[i];
      pth->ClrAdr = &tab[ i * par->ClrLinSiz ];

      pthread_mutex_lock(&pth->mtx);
      pthread_cond_signal(&pth->cnd);
      pthread_mutex_unlock(&pth->mtx);
   }

   // Wait for each thread to complete
   while(par->WrkCpt < par->NmbCpu)
      pthread_cond_wait(&par->ParCnd, &par->ParMtx);

   pthread_mutex_unlock(&par->ParMtx);

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Wait for a condition, launch and detach a user procedure                   */
/*----------------------------------------------------------------------------*/

int LaunchPipeline(  int64_t ParIdx, void *prc,
                     void *PtrArg, int NmbDep, int *DepTab )
{
   int i;
   PipSct *NewPip=NULL;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance and the number of pipes and dependencies
   if( !ParIdx || (NmbDep > MaxPipDep) || (par->NmbPip >= MaxTotPip) )
      return(0);

   // Allocate and setup a new pipe
   if(!(NewPip = calloc(1, sizeof(PipSct))))
      return(0);

   NewPip->prc = prc;
   NewPip->arg = PtrArg;
   NewPip->par = par;
   NewPip->NmbDep = NmbDep;

   for(i=0;i<NmbDep;i++)
      NewPip->DepTab[i] = DepTab[i];

   // In case variable arguments where passed through the common LPlib structure
   // Copy them in the pipeline's own arguments table
   if(par->NmbF77Arg)
   {
      NewPip->NmbF77Arg = par->NmbF77Arg;

      for(i=0;i<par->NmbF77Arg;i++)
         NewPip->F77ArgTab[i] = par->F77ArgTab[i];
   }

   // Lock pipe mutex, increment pipe counter
   // and launch the pipe regardless dependencies
   pthread_mutex_lock(&par->PipMtx);
   NewPip->idx = ++par->NmbPip;
   par->PenPip++;

   if(par->StkSiz)
   {
      pthread_attr_init(&NewPip->atr);
      NewPip->StkSiz = par->StkSiz;
      NewPip->UsrStk = malloc(par->StkSiz);
      pthread_attr_setstackaddr(&NewPip->atr, NewPip->UsrStk);
      pthread_attr_setstacksize(&NewPip->atr, NewPip->StkSiz);
      //pthread_attr_setstack(&NewPip->atr, NewPip->UsrStk, NewPip->StkSiz);
      pthread_create(&NewPip->pth, &NewPip->atr, PipHdl, (void *)NewPip);
   }
   else
      pthread_create(&NewPip->pth, NULL, PipHdl, (void *)NewPip);

   // Make the thread unjoinable to work around a memory leak in some systems
   pthread_detach(NewPip->pth);
   pthread_mutex_unlock(&par->PipMtx);

   return(NewPip->idx);
}


/*----------------------------------------------------------------------------*/
/* Launch a pipeline procedure with variable numbers of arguments             */
/*----------------------------------------------------------------------------*/

int LaunchPipelineMultiArg(int64_t ParIdx, int NmbDep, int *DepTab, 
                           void *prc, int NmbArg, ... )
{
   int i, ret;
   va_list ArgLst;
   ParSct *par = (ParSct *)ParIdx;

   if(NmbArg > MaxF77Arg)
      return(-1);

   // Get the variable arguments and store then in the common LPlib structure
   // This is a temporary storage as the will be copied to the pipeline's own
   // arguments table after it has been allocated by LaunchPipeline()
   par->NmbF77Arg = NmbArg;
   va_start(ArgLst, NmbArg);

   for(i=0;i<NmbArg;i++)
      par->F77ArgTab[i] = va_arg(ArgLst, void *);

   va_end(ArgLst);

   ret = LaunchPipeline(ParIdx, prc, NULL, NmbDep, DepTab);

   // Clear the temporary arguments table
   par->NmbF77Arg = 0;

   return(ret);
}


/*----------------------------------------------------------------------------*/
/* Thread handler launching and waitint for user's procedure                  */
/*----------------------------------------------------------------------------*/

static void *PipHdl(void *ptr)
{
   int RunFlg=0, i;
   PipSct *pip = (PipSct *)ptr;
   ParSct *par = pip->par;
   void (*prc)(void *);

   // Wait for conditions to be met
   do
   {
      pthread_mutex_lock(&par->PipMtx);

      if(par->RunPip < par->NmbCpu)
      {
         RunFlg = 1;

         for(i=0;i<pip->NmbDep;i++)
            if(!GetBit(par->PipWrd, pip->DepTab[i]))
            {
               RunFlg = 0;
               break;
            }
      }

      if(!RunFlg)
      {
         pthread_mutex_unlock(&par->PipMtx);
#ifdef _WIN32
         Sleep(1);
#else
         usleep(1000);
#endif
      }
   }while(!RunFlg);

   // Execute the user's procedure and set the flag to 2 (done)
   prc = (void (*)(void *))pip->prc;
   par->RunPip++;

   pthread_mutex_unlock(&par->PipMtx);

   if(pip->NmbF77Arg)
      CalF77Pip(pip, pip->prc);
   else
      prc(pip->arg);

   pthread_mutex_lock(&par->PipMtx);
   SetBit(par->PipWrd, pip->idx);
   par->PenPip--;
   par->RunPip--;
   free(pip);
   pthread_mutex_unlock(&par->PipMtx);

   return(NULL);
}


/*----------------------------------------------------------------------------*/
/* Wait for all pipelined procedures to complete                              */
/*----------------------------------------------------------------------------*/

void WaitPipeline(int64_t ParIdx)
{
   int PenPip;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   do
   {
      pthread_mutex_lock(&par->PipMtx);
      PenPip = par->PenPip;
      pthread_mutex_unlock(&par->PipMtx);
#ifdef _WIN32
      Sleep(1);
#else
      usleep(1000);
#endif
   }while(PenPip);
}


/*----------------------------------------------------------------------------*/
/* Wait for all pipelined procedures to complete                              */
/*----------------------------------------------------------------------------*/

void GetLplibInformation(int64_t ParIdx, int *NmbCpu, int *NmbTyp)
{
   int i;
   ParSct *par = (ParSct *)ParIdx;

   // Get and check lib parallel instance
   if(!ParIdx)
      return;

   *NmbCpu = par->NmbCpu;
   *NmbTyp = 0;

   // Count the number of used types
   for(i=1;i<=MaxTyp;i++)
      if(par->TypTab[i].NmbLin)
         (*NmbTyp)++;
}


/*----------------------------------------------------------------------------*/
/* Call qsort in parallel then call merge sort in serial                      */
/*----------------------------------------------------------------------------*/

void ParallelQsort(  int64_t ParIdx, void *base, size_t nel, size_t width,
                     int (*compar)(const void *, const void *) )
{
   (void)(ParIdx);
   qsort(base, nel, width, compar);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates                               */
/*----------------------------------------------------------------------------*/

static void RenPrc(itg BegIdx, itg EndIdx, int PthIdx, ArgSct *arg)
{
   uint64_t IntCrd[3], m=1ULL<<63, cod;
   itg i;
   int j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
   double dbl;
   int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  // Z curve = {5,4,7,6,1,0,3,2}
   int HilCod[8][8] = {
      {0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1},
      {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},
      {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7},
      {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7} };
   (void)(PthIdx);

   for(i=BegIdx; i<=EndIdx; i++)
   {
      // Convert double precision coordinates to integers
      for(j=0;j<3;j++)
      {
         dbl = (arg->crd[3*i + j] - arg->box[j]) * arg->box[j+3];
         IntCrd[j] = (uint64_t)dbl;
      }

      // Binary hilbert renumbering loop
      cod = 0;

      for(j=0;j<8;j++)
         rot[j] = GeoCod[j];

      for(b=0;b<21;b++)
      {
         GeoWrd = 0;

         for(j=0;j<3;j++)
         {
            if(IntCrd[j] & m)
               GeoWrd |= BitTab[j];

            IntCrd[j] = IntCrd[j]<<1;
         }

         NewWrd = rot[ GeoWrd ];
         cod = cod<<3 | NewWrd;

         for(j=0;j<8;j++)
            rot[j] = HilCod[ NewWrd ][ rot[j] ];
      }
      arg->idx[2*i+0] = cod;
      arg->idx[2*i+1] = i;
   }
}


/*----------------------------------------------------------------------------*/
/* Comparison of two items for the qsort                                      */
/*----------------------------------------------------------------------------*/

int CmpPrc(const void *a, const void *b)
{
   uint64_t *pa = (uint64_t *)a, *pb = (uint64_t *)b;
   return(*pa > *pb ? 1 : -1);
}

void PipSrt(PipArgSct *arg)
{
   qsort(arg->base, arg->nel, arg->width, arg->compar);
}


/*----------------------------------------------------------------------------*/
/* Renumber a set of coordinates through a Hilbert SFC                        */
/*----------------------------------------------------------------------------*/

int HilbertRenumbering( int64_t ParIdx, itg NmbLin, double box[6],
                        double *crd, uint64_t *idx )
{
   itg i;
   int j, NewTyp, stat[ (1<<HshBit)+1 ], NmbPip;
   uint64_t bound[ MaxPth ][2], cpt, sum, (*tab)[2];
   size_t NmbByt;
   double len = pow(2,64);
   ParSct *par = (ParSct *)ParIdx;
   ArgSct arg;
   PipArgSct PipArg[ MaxPth ];

   // Get and check lib parallel instance
   if(!ParIdx)
     return(0);

   // Setup the bounding box and a data type,
   // then give a Hilbert code to each entries
   NewTyp = NewType(ParIdx, NmbLin);
   arg.crd = crd;
   arg.idx = idx;
   arg.box[0] = box[0];
   arg.box[1] = box[1];
   arg.box[2] = box[2];
   arg.box[3] = len / (box[3] - box[0]);
   arg.box[4] = len / (box[4] - box[1]);
   arg.box[5] = len / (box[5] - box[2]);


   if(NmbLin < 10000)
   {
      RenPrc(1, NmbLin, 0, (void *)&arg);

      qsort(&idx[2*1 + 0], NmbLin, 2 * sizeof(int64_t), CmpPrc);

      for(i=1;i<=NmbLin;i++)
         idx[ 2*idx[2*i + 1] + 0] = i;
   }
   else
   {
      LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc, (void *)&arg);

      for(i=0;i<1<<HshBit;i++)
         stat[i] = 0;

      for(i=1;i<=NmbLin;i++)
         stat[ idx[2*i + 0] >> (64 - HshBit) ]++;

      for(i=0;i<MaxPth;i++)
         bound[i][0] = bound[i][1] = 0;

      NmbPip = 0;
      sum = cpt = 0;
   
      for(i=0;i<1<<HshBit;i++)
      {
         bound[ NmbPip ][0] += stat[i];
         cpt++;

         if(bound[ NmbPip ][0] >= (size_t)NmbLin / par->NmbCpu)
         {
            bound[ NmbPip ][1] = cpt << (64 - HshBit);
            NmbPip++;
         }
      }

      bound[ NmbPip ][1] = 1LL<<63;
   
      NmbByt = NmbLin+1;
      NmbByt *= 2 * sizeof(int64_t);
   
      tab = malloc(NmbByt);

      if(!tab)
         return(0);

      for(i=0;i<=NmbPip;i++)
      {
         cpt = bound[i][0];
         bound[i][0] = sum;
         sum += cpt;

         PipArg[i].base = &tab[ bound[i][0] ];
         PipArg[i].nel = cpt;
         PipArg[i].width = 2 * sizeof(int64_t);
         PipArg[i].compar = CmpPrc;
      }

      for(i=1;i<=NmbLin;i++)
         for(j=0;j<=NmbPip;j++)
            if(idx[2*i + 0] <= bound[j][1])
            {
               tab[ bound[j][0] ][0] = idx[2*i + 0];
               tab[ bound[j][0] ][1] = i;
               bound[j][0]++;
               break;
            }

      for(i=0;i<NmbPip;i++)
         LaunchPipeline(ParIdx, PipSrt, &PipArg[i], 0, NULL);

      WaitPipeline(ParIdx);

      for(i=1;i<=NmbLin;i++)
      {
         idx[2*i + 1] = tab[i-1][1];
         idx[ 2*tab[i-1][1] + 0] = i;
      }

      free(tab);
   }

   return(1);
}


/*----------------------------------------------------------------------------*/
/* Compute the hilbert code from 2d coordinates                               */
/*----------------------------------------------------------------------------*/

static void RenPrc2D(itg BegIdx, itg EndIdx, int PthIdx, ArgSct *arg)
{
   uint64_t IntCrd[2], m=1ULL<<62, cod;
   itg i;
   int j, b, GeoWrd, NewWrd, BitTab[2] = {1,2};
   double dbl;
   int rot[4], GeoCod[4]={1,2,0,3};
   int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
   (void)(PthIdx);

   for(i=BegIdx; i<=EndIdx; i++)
   {
      // Convert double precision coordinates to integers
      for(j=0;j<2;j++)
      {
         dbl = (arg->crd2[3*i + j] - arg->box[j]) * arg->box[j+2];
         IntCrd[j] = (uint64_t)dbl;
      }

      // Binary hilbert renumbering loop
      cod = 0;

      for(j=0;j<4;j++)
         rot[j] = GeoCod[j];

      for(b=0;b<31;b++)
      {
         GeoWrd = 0;

         for(j=0;j<2;j++)
         {
            if(IntCrd[j] & m)
               GeoWrd |= BitTab[j];

            IntCrd[j] = IntCrd[j]<<1;
         }

         NewWrd = rot[ GeoWrd ];
         cod = cod<<2 | NewWrd;

         for(j=0;j<4;j++)
            rot[j] = HilCod[ NewWrd ][ rot[j] ];
      }

      arg->idx[2*i + 0] = cod;
      arg->idx[2*i + 1] = i;
   }
}


/*----------------------------------------------------------------------------*/
/* Renumber a set of 2D coordinates through a Hilbert SFC                     */
/*----------------------------------------------------------------------------*/

//int HilbertRenumbering2D(  int64_t ParIdx, itg NmbLin, double box[4],
//                           double (*crd)[2], uint64_t (*idx)[2] )
//{
//   itg i;
//   int NewTyp;
//   double len = pow(2,62);
//   ArgSct arg;
//
//   // Get and check lib parallel instance
//   if(!ParIdx)
//      return(0);
//
//   NewTyp = NewType(ParIdx, NmbLin);
//   arg.crd2 = crd;
//   arg.idx = idx;
//   arg.box[0] = box[0];
//   arg.box[1] = box[1];
//   arg.box[2] = len / (box[2] - box[0]);
//   arg.box[3] = len / (box[3] - box[1]);
//
//   LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc2D, (void *)&arg);
//   ParallelQsort(ParIdx, &idx[1][0], NmbLin, 2 * sizeof(int64_t), CmpPrc);
//
//   for(i=1;i<=NmbLin;i++)
//      idx[ idx[i][1] ][0] = i;
//
//   return(1);
//}


/*----------------------------------------------------------------------------*/
/* Starts or stops the given timer                                            */
/*----------------------------------------------------------------------------*/

double GetWallClock()
{
#ifdef _WIN32
   struct __timeb64 tb;
   _ftime64(&tb);
   return((double)tb.time + (double)tb.millitm/1000.);
#else
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return(tp.tv_sec + tp.tv_usec / 1000000.);
#endif
}


/*----------------------------------------------------------------------------*/
/* Duplication macros                                                        */
/*----------------------------------------------------------------------------*/

#define DUP(s,n) DUP ## n (s)
#define DUP1(s) s
#define DUP2(s) DUP1(s),s
#define DUP3(s) DUP2(s),s
#define DUP4(s) DUP3(s),s
#define DUP5(s) DUP4(s),s
#define DUP6(s) DUP5(s),s
#define DUP7(s) DUP6(s),s
#define DUP8(s) DUP7(s),s
#define DUP9(s) DUP8(s),s
#define DUP10(s) DUP9(s),s
#define DUP11(s) DUP10(s),s
#define DUP12(s) DUP11(s),s
#define DUP13(s) DUP12(s),s
#define DUP14(s) DUP13(s),s
#define DUP15(s) DUP14(s),s
#define DUP16(s) DUP15(s),s
#define DUP17(s) DUP16(s),s
#define DUP18(s) DUP17(s),s
#define DUP19(s) DUP18(s),s
#define DUP20(s) DUP19(s),s


#define ARG(a,n) ARG ## n (a)
#define ARG1(a) a[0]
#define ARG2(a) ARG1(a),a[1]
#define ARG3(a) ARG2(a),a[2]
#define ARG4(a) ARG3(a),a[3]
#define ARG5(a) ARG4(a),a[4]
#define ARG6(a) ARG5(a),a[5]
#define ARG7(a) ARG6(a),a[6]
#define ARG8(a) ARG7(a),a[7]
#define ARG9(a) ARG8(a),a[8]
#define ARG10(a) ARG9(a),a[9]
#define ARG11(a) ARG10(a),a[10]
#define ARG12(a) ARG11(a),a[11]
#define ARG13(a) ARG12(a),a[12]
#define ARG14(a) ARG13(a),a[13]
#define ARG15(a) ARG14(a),a[14]
#define ARG16(a) ARG15(a),a[15]
#define ARG17(a) ARG16(a),a[16]
#define ARG18(a) ARG17(a),a[17]
#define ARG19(a) ARG18(a),a[18]
#define ARG20(a) ARG19(a),a[19]


/*----------------------------------------------------------------------------*/
/* Call a fortran thread with 1 to 20 arguments                               */
/*----------------------------------------------------------------------------*/

static void CalF77Prc(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   switch(par->NmbF77Arg)
   {
      case 1 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 1)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 1)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 1));
      }break;

      case 2 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 2)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 2)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 2));
      }break;

      case 3 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 3)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 3)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 3));
      }break;

      case 4 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 4)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 4)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 4));
      }break;

      case 5 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 5)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 5)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 5));
      }break;

      case 6 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 6)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 6)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 6));
      }break;

      case 7 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 7)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 7)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 7));
      }break;

      case 8 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 8)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 8)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 8));
      }break;

      case 9 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 9)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 9)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 9));
      }break;

      case 10 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 10)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 10)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 10));
      }break;

      case 11 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 11)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 11)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 11));
      }break;

      case 12 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 12)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 12)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 12));
      }break;

      case 13 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 13)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 13)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 13));
      }break;

      case 14 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 14)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 14)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 14));
      }break;

      case 15 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 15)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 15)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 15));
      }break;

      case 16 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 16)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 16)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 16));
      }break;

      case 17 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 17)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 17)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 17));
      }break;

      case 18 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 18)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 18)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 18));
      }break;

      case 19 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 19)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 19)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 19));
      }break;

      case 20 :
      {
         void (*prc1)(itg *, itg *, int *, DUP(void *, 20)) =
            (void (*)(itg *, itg *, int *, DUP(void *, 20)))par->prc;
         prc1(&BegIdx, &EndIdx, &PthIdx, ARG(par->F77ArgTab, 20));
      }break;
   }
}


/*----------------------------------------------------------------------------*/
/* Call a C thread with 1 to 20 arguments                                     */
/*----------------------------------------------------------------------------*/

static void CalVarArgPrc(itg BegIdx, itg EndIdx, int PthIdx, ParSct *par)
{
   switch(par->NmbVarArg)
   {
      case 1 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 1)) =
            (void (*)(itg, itg, int, DUP(void *, 1)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 1));
      }break;

      case 2 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 2)) =
            (void (*)(itg, itg, int, DUP(void *, 2)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 2));
      }break;

      case 3 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 3)) =
            (void (*)(itg, itg, int, DUP(void *, 3)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 3));
      }break;

      case 4 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 4)) =
            (void (*)(itg, itg, int, DUP(void *, 4)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 4));
      }break;

      case 5 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 5)) =
            (void (*)(itg, itg, int, DUP(void *, 5)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 5));
      }break;

      case 6 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 6)) =
            (void (*)(itg, itg, int, DUP(void *, 6)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 6));
      }break;

      case 7 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 7)) =
            (void (*)(itg, itg, int, DUP(void *, 7)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 7));
      }break;

      case 8 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 8)) =
            (void (*)(itg, itg, int, DUP(void *, 8)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 8));
      }break;

      case 9 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 9)) =
            (void (*)(itg, itg, int, DUP(void *, 9)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 9));
      }break;

      case 10 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 10)) =
            (void (*)(itg, itg, int, DUP(void *, 10)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 10));
      }break;

      case 11 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 11)) =
            (void (*)(itg, itg, int, DUP(void *, 11)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 11));
      }break;

      case 12 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 12)) =
            (void (*)(itg, itg, int, DUP(void *, 12)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 12));
      }break;

      case 13 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 13)) =
            (void (*)(itg, itg, int, DUP(void *, 13)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 13));
      }break;

      case 14 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 14)) =
            (void (*)(itg, itg, int, DUP(void *, 14)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 14));
      }break;

      case 15 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 15)) =
            (void (*)(itg, itg, int, DUP(void *, 15)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 15));
      }break;

      case 16 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 16)) =
            (void (*)(itg, itg, int, DUP(void *, 16)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 16));
      }break;

      case 17 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 17)) =
            (void (*)(itg, itg, int, DUP(void *, 17)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 17));
      }break;

      case 18 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 18)) =
            (void (*)(itg, itg, int, DUP(void *, 18)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 18));
      }break;

      case 19 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 19)) =
            (void (*)(itg, itg, int, DUP(void *, 19)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 19));
      }break;

      case 20 :
      {
         void (*prc1)(itg, itg, int, DUP(void *, 20)) =
            (void (*)(itg, itg, int, DUP(void *, 20)))par->prc;
         prc1(BegIdx, EndIdx, PthIdx, ARG(par->F77ArgTab, 20));
      }break;
   }
}


/*----------------------------------------------------------------------------*/
/* Call a fortran pipeline with 1 to 20 arguments                             */
/*----------------------------------------------------------------------------*/

static void CalF77Pip(PipSct *pip, void *prc)
{
   switch(pip->NmbF77Arg)
   {
      case 1 :
      {
         void (*prc1)(DUP(void *, 1)) = (void (*)(DUP(void *, 1)))prc;
         prc1(ARG(pip->F77ArgTab, 1));
      }break;

      case 2 :
      {
         void (*prc1)(DUP(void *, 2)) = (void (*)(DUP(void *, 2)))prc;
         prc1(ARG(pip->F77ArgTab, 2));
      }break;

      case 3 :
      {
         void (*prc1)(DUP(void *, 3)) = (void (*)(DUP(void *, 3)))prc;
         prc1(ARG(pip->F77ArgTab, 3));
      }break;

      case 4 :
      {
         void (*prc1)(DUP(void *, 4)) = (void (*)(DUP(void *, 4)))prc;
         prc1(ARG(pip->F77ArgTab, 4));
      }break;

      case 5 :
      {
         void (*prc1)(DUP(void *, 5)) = (void (*)(DUP(void *, 5)))prc;
         prc1(ARG(pip->F77ArgTab, 5));
      }break;

      case 6 :
      {
         void (*prc1)(DUP(void *, 6)) = (void (*)(DUP(void *, 6)))prc;
         prc1(ARG(pip->F77ArgTab, 6));
      }break;

      case 7 :
      {
         void (*prc1)(DUP(void *, 7)) = (void (*)(DUP(void *, 7)))prc;
         prc1(ARG(pip->F77ArgTab, 7));
      }break;

      case 8 :
      {
         void (*prc1)(DUP(void *, 8)) = (void (*)(DUP(void *, 8)))prc;
         prc1(ARG(pip->F77ArgTab, 8));
      }break;

      case 9 :
      {
         void (*prc1)(DUP(void *, 9)) = (void (*)(DUP(void *, 9)))prc;
         prc1(ARG(pip->F77ArgTab, 9));
      }break;

      case 10 :
      {
         void (*prc1)(DUP(void *, 10)) = (void (*)(DUP(void *, 10)))prc;
         prc1(ARG(pip->F77ArgTab, 10));
      }break;

      case 11 :
      {
         void (*prc1)(DUP(void *, 11)) = (void (*)(DUP(void *, 11)))prc;
         prc1(ARG(pip->F77ArgTab, 11));
      }break;

      case 12 :
      {
         void (*prc1)(DUP(void *, 12)) = (void (*)(DUP(void *, 12)))prc;
         prc1(ARG(pip->F77ArgTab, 12));
      }break;

      case 13 :
      {
         void (*prc1)(DUP(void *, 13)) = (void (*)(DUP(void *, 13)))prc;
         prc1(ARG(pip->F77ArgTab, 13));
      }break;

      case 14 :
      {
         void (*prc1)(DUP(void *, 14)) = (void (*)(DUP(void *, 14)))prc;
         prc1(ARG(pip->F77ArgTab, 14));
      }break;

      case 15 :
      {
         void (*prc1)(DUP(void *, 15)) = (void (*)(DUP(void *, 15)))prc;
         prc1(ARG(pip->F77ArgTab, 15));
      }break;

      case 16 :
      {
         void (*prc1)(DUP(void *, 16)) = (void (*)(DUP(void *, 16)))prc;
         prc1(ARG(pip->F77ArgTab, 16));
      }break;

      case 17 :
      {
         void (*prc1)(DUP(void *, 17)) = (void (*)(DUP(void *, 17)))prc;
         prc1(ARG(pip->F77ArgTab, 17));
      }break;

      case 18 :
      {
         void (*prc1)(DUP(void *, 18)) = (void (*)(DUP(void *, 18)))prc;
         prc1(ARG(pip->F77ArgTab, 18));
      }break;

      case 19 :
      {
         void (*prc1)(DUP(void *, 19)) = (void (*)(DUP(void *, 19)))prc;
         prc1(ARG(pip->F77ArgTab, 19));
      }break;

      case 20 :
      {
         void (*prc1)(DUP(void *, 20)) = (void (*)(DUP(void *, 20)))prc;
         prc1(ARG(pip->F77ArgTab, 20));
      }break;
   }
}


/*----------------------------------------------------------------------------*/
/* Fortran 77 API                                                             */
/*----------------------------------------------------------------------------*/

/*
#ifndef __cplusplus
   #if defined(F77_NO_UNDER_SCORE)
   #define call(x) x
   #else
   #define call(x) x ## _
   #endif
   
   
   int64_t call(initparallel)(int *NmbCpu)
   {
      return(InitParallel(*NmbCpu));
   }
   
   void call(stopparallel)(int64_t *ParIdx)
   {
      StopParallel(*ParIdx);
   }
   
   int call(newtype)(int64_t *ParIdx, itg *NmbLin)
   {
      return(NewType(*ParIdx, *NmbLin));
   }
   
   void call(freetype)(int64_t *ParIdx, int *TypIdx)
   {
      FreeType(*ParIdx, *TypIdx);
   }
   
   int call(begindependency)(int64_t *ParIdx, int *TypIdx1, int *TypIdx2)
   {
      return(BeginDependency(*ParIdx, *TypIdx1, *TypIdx2));
   }
   
   
   void call(adddependency)(int64_t *ParIdx, itg *idx1, itg *idx2)
   {
      AddDependency(*ParIdx, *idx1, *idx2);
   }
   
   
   void call(enddependency)(int64_t *ParIdx, float *DepSta)
   {
      EndDependency(*ParIdx, DepSta);
   }
   
   
   float call( launchparallel)(int64_t *ParIdx, int *typ1,
               int *typ2, void *prc, int *NmbArg, ... )
   {
      int i;
      float acc;
      va_list ArgLst;
      ParSct *par = (ParSct *)*ParIdx;
   
      par->NmbF77Arg = *NmbArg;
      va_start(ArgLst, NmbArg);
   
      for(i=0;i<*NmbArg;i++)
         par->F77ArgTab[i] = va_arg(ArgLst, void *);
   
      va_end(ArgLst);
   
      acc = LaunchParallel(*ParIdx, *typ1, *typ2, prc, NULL);
      par->NmbF77Arg = 0;
   
      return(acc);
   }
   
   
   void call(parallelmemclear)(int64_t *ParIdx, void *ptr, size_t *siz)
   {
      ParallelMemClear(*ParIdx, ptr, (size_t)*siz);
   }
   
   
   void call(getlplibinformation)(int64_t *ParIdx, int *NmbCpu, int *NmbTyp)
   {
      GetLplibInformation(*ParIdx, NmbCpu, NmbTyp);
   }
   
   int call(launchpipeline)(  int64_t *ParIdx, int *NmbDep,
                              int *DepTab, void *prc, int *NmbArg, ... )
   {
      int i, ret;
      va_list ArgLst;
      ParSct *par = (ParSct *)*ParIdx;
   
      par->NmbF77Arg = *NmbArg;
      va_start(ArgLst, NmbArg);
   
      for(i=0;i<*NmbArg;i++)
         par->F77ArgTab[i] = va_arg(ArgLst, void *);
   
      va_end(ArgLst);
      ret = LaunchPipeline(*ParIdx, prc, NULL, *NmbDep, DepTab);
      return(ret);
   }
   
   void call(waitpipeline)(int64_t *ParIdx)
   {
      WaitPipeline(*ParIdx);
   }
   
   double call(hilbertrenumbering)( int64_t *ParIdx, itg *NmbLin, double box[6],
                                    double (*crd)[3], uint64_t *Old2New)
   {
      int i;
      uint64_t (*idx)[2];
   
      // Alloc a temporary renumbering index table
      if(!(idx = malloc((size_t)(*NmbLin+1) * 2 * sizeof(int64_t))))
        return(0.);
   
      // First entry of crd should not be used, crd(1) in F77 is crd[0] in C
      if(!HilbertRenumbering((int)*ParIdx, (int)*NmbLin, box, &crd[-1], idx))
        return(0.);
   
      // Copy the old to new only in the fortran table and free the temp idx table
      for(i=1;i<=*NmbLin;i++)
         Old2New[i-1] = idx[i][0];
   
      free(idx);
   
      return(1.);
   }
   
   void call(hilbertrenumbering2d)( int64_t *ParIdx, itg *NmbLin, double box[4],
                                    double (*crd)[2], uint64_t (*idx)[2])
   {
      HilbertRenumbering2D(*ParIdx, *NmbLin, box, &crd[-1], &idx[-1]);
   }
   
   int call(getnumberofcores)()
   {
      return(GetNumberOfCores());
   }
   
   double call(getwallclock)()
   {
      return(GetWallClock());
   }
#endif
*/
