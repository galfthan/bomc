/*
BOMC - Monte Carlo based generator of defect-free amorphous samples of Si and SiO2

If you publish research done using Bomc then you are kindly asked to cite: 
 S. von Alfthan, A. Kuronen, and K. Kaski, Phys. Rev. B 68, 073203 (2003).


Copyright Sebastian von Alfthan  galfthan@iki.fi 2002, 2003, 2014.


 This file is part of Bomc.	 

Bomc is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.  

Bomc is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Bomc.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include "shared.h"
#include "conjgrad.h"
#include "fileio.h"
#include "nrMinSubs.h"
#include "parallel.h"

#define ITMAX 10000
#define EPS 1.0e-10
/*conjgrad based on NR version */

struct optimPar{
  int useOptim;
  int startIteration;
  double c;
  int checkInterval;
};

struct fullCgPar{
  struct optimPar optim;
  int allowAdaptive;
  double tolerance;
};
struct partCgPar{
  struct optimPar optim;
  int allowAdaptive;
  double tolerance;
  int enlargeInterval;
  int enlargements;
  double enlTreshold;

};

struct initPartPar{
  int allowAdaptive;
  struct partCgPar part;
  int doInitPart;
  int doPartCg;
  int doFullCg;
};
 
struct cgPar{
  struct fullCgPar full;
  struct partCgPar part;
  struct initPartPar initPart;
};

static struct cgPar cgPar;
static  double meanPartSize=-1;
static double meanInitPartSize=-1;

void  getCgStatus(char *message){
  sprintf(message,"Average part size: %.1f initPartSize: %.1f", meanPartSize,meanInitPartSize);
}


void initCgPar(void){
  static int firsttime=1;
  int error=0;
  
  if(firsttime)
    firsttime=0;
  else
    return;
  
  
  openParameterFile("cg.par");
  
  /*full cg */
  error+=getParValue("full_tolerance"          ,&cgPar.full.tolerance          ,"%lf");
  error+=getParValue("full_allowAdaptive"      ,&cgPar.full.allowAdaptive      ,"%d");

  error+=getParValue("full_optim_on"            ,&cgPar.full.optim.useOptim      ,"%d");
  error+=getParValue("full_optim_startIteration",&cgPar.full.optim.startIteration,"%d");
  error+=getParValue("full_optim_c"             ,&cgPar.full.optim.c             ,"%lf");
  error+=getParValue("full_optim_checkInterval" ,&cgPar.full.optim.checkInterval ,"%d");
  
  /*fpart cg */
  error+=getParValue("part_tolerance"          ,&cgPar.part.tolerance          ,"%lf");
  error+=getParValue("part_allowAdaptive"      ,&cgPar.part.allowAdaptive      ,"%d");
  error+=getParValue("part_enlargeInterval"    ,&cgPar.part.enlargeInterval    ,"%d");
  error+=getParValue("part_enlargements"       ,&cgPar.part.enlargements       ,"%d");
  error+=getParValue("part_enlTreshold"        ,&cgPar.part.enlTreshold        ,"%lf");
  
  error+=getParValue("part_optim_on"           ,&cgPar.part.optim.useOptim     ,"%d");
  error+=getParValue("part_optim_c"            ,&cgPar.part.optim.c            ,"%lf");
  error+=getParValue("part_optim_checkInterval",&cgPar.part.optim.checkInterval,"%d");
  
  /*initial optimization in part cg */
  error+=getParValue("initPart_doInitPart"             ,&cgPar.initPart.doInitPart               ,"%d");
  error+=getParValue("initPart_doPartCg"                ,&cgPar.initPart.doPartCg                ,"%d");
  error+=getParValue("initPart_doFullCg"                ,&cgPar.initPart.doFullCg                ,"%d");

  error+=getParValue("initpartpart_allowAdaptive"      ,&cgPar.initPart.allowAdaptive            ,"%d");  
  error+=getParValue("initPart_part_enlargeInterval"    ,&cgPar.initPart.part.enlargeInterval    ,"%d");
  error+=getParValue("initPart_part_enlargements"       ,&cgPar.initPart.part.enlargements       ,"%d");
  error+=getParValue("initPart_part_enlTreshold"        ,&cgPar.initPart.part.enlTreshold        ,"%lf");
  
  error+=getParValue("initPart_part_optim_on"           ,&cgPar.initPart.part.optim.useOptim     ,"%d");
  error+=getParValue("initPart_part_optim_c"            ,&cgPar.initPart.part.optim.c            ,"%lf");
  error+=getParValue("initPart_part_optim_checkInterval",&cgPar.initPart.part.optim.checkInterval,"%d");
  
  closeParameterFile();
  
  cgPar.initPart.part.optim.startIteration=0; /*the optimization in the part cg start when part is at final size... */
  cgPar.part.optim.startIteration=0; /*the optimization in the part cg start when part is at final size... */
  cgPar.initPart.part.tolerance=0.0; /*this is not used, the intial part is always the same length */
  
  cgPar.initPart.part.enlTreshold=POW2(cgPar.initPart.part.enlTreshold);
  cgPar.part.enlTreshold=POW2(cgPar.part.enlTreshold);
 
  
  if(error){
    if(getRank()==0) printf("could not read needed parameters for CG from cgnew.par\n");
    exit(0);
  }
  if(cgPar.initPart.doPartCg==0 && cgPar.initPart.doFullCg==0){
    if(getRank()==0) printf("Either initPart_doPartCg or initPart_doFullCg have to be on\n");
    exit(0);
  }
		
}



double fullCg(struct systemPos pos,struct parameters *par,int allowOptim,int isHarmonic,double accPotVal,double (*pot)(struct parameters*,struct systemPos,int calcForces)){
  
  int j,its;
  double gg,gam,fp,dgg;
  double totForce;
  double potVal;
  
  static double *g,*h,*xi;
  static int firsttime=1;
  static int numOfCalls=0;
  
  int tryOptimization;

  int checkOptim;
  double checkOptimPotential;
  double checkOptimForce;
  int checkOptimAborted=0;
  int parallelAbortion=0;
  
  if(firsttime){
    firsttime=0;
	 
    g=malloc(sizeof(double)*pos.nAtoms*3);
    h=malloc(sizeof(double)*pos.nAtoms*3);
    if(g==NULL|| h==NULL){
      printf("out of memory\n");
      exit(2);
    }
    initCgPar();
  }
  
  
  numOfCalls++;
  checkOptim=(numOfCalls%cgPar.full.optim.checkInterval)==0;
  
  
  fp=(*pot)(par,pos,1);
  potVal=fp;
  xi=pos.xa;
  
  for (j=0;j<3*pos.nAtoms;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  
  for (its=0;its<ITMAX;its++) {
    linmin(&pos,h,par,&potVal,cgPar.full.allowAdaptive,pot);
	 if (  2.0*fabs(potVal-fp) <= cgPar.full.tolerance*(fabs(potVal)+fabs(fp)+EPS)) 
      break;
    
    fp=(*pot)(par,pos,1); /*xi=pos.xa updated*/
    
    dgg=gg=totForce=0.0;
    for (j=0;j<3*pos.nAtoms;j++){
      totForce+=xi[j]*xi[j];
      gg      += g[j]*g[j];
      dgg     += (xi[j]+g[j])*xi[j];
    }
    
    if (gg == 0.0) 
      break;
    
    tryOptimization=allowOptim && cgPar.full.optim.useOptim && (its>=cgPar.full.optim.startIteration) ;
	 
    if(tryOptimization  &&  fp-cgPar.full.optim.c*totForce>accPotVal){
      if(checkOptim){
		  if(!checkOptimAborted){   /*onlu tke first one */
			 checkOptimAborted=1;		  /* ok tried to break, used later on to see if c value was alright or not */
			 checkOptimPotential=fp;
			 checkOptimForce=totForce;
		  }
      }
      else{
		  if(DEBUG) printf("(%d) fullCg aborted on step %d\n",getRank(),its);
		  break;
      }
    }
#ifdef PARALLEL
    { /*this block aborts this minimization if some oter process has an accepted step */
      int flag;
      MPI_Status status;
      MPI_Iprobe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&flag,&status); /*probe if somebody has sent a message, i.e. if a step has been accepted*/
      if(flag){
		  parallelAbortion=1;
		  potVal=1e30; /* do not accept this step*/
		  break;
      }
    }
#endif
    
    gam=dgg/gg;
    for (j=0;j<3*pos.nAtoms;j++){
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  
  /*check optimization*/
  if(!parallelAbortion && allowOptim && cgPar.full.optim.useOptim && checkOptim && checkOptimAborted){
    double cGuess=1.1*(checkOptimPotential-potVal)/checkOptimForce; /*corrects c value to a more realistic one.. */
    if(cgPar.full.optim.c<cGuess)
      cgPar.full.optim.c=cGuess;
    
    if(DEBUG) printf("(%d) cgpar_full_optim_c:%g  \n",getRank(),cgPar.full.optim.c);
  }
  
  if(DEBUG) printf("(%d) %d steps in fullCg, pot=%g\n",getRank(),its,potVal/pos.nAtoms*TOEV);
  return potVal;
}  


double partCg(struct systemPos pos,struct parameters *par,int *part, int *atomsInPart,double accPotVal,	
	      double (*partpot)(struct parameters*,struct systemPos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold),
	      double (*pot)(struct parameters*,struct systemPos,int calcForces))
{
  int j,jps,d,its;
  double gg,gam,fp,dgg,totForce;
 
  double potVal;
  double restPot;
  
  static double *g,*h,*xi;
  static int firsttime=1;
  
  int restPotCalculated=0;
  int partAtFinalSize;
  
  static int numOfCalls=0;

  int checkOptim;
  double checkOptimPotential;
  double checkOptimForce;
  int checkOptimAborted=0;
  int parallelAbortion=0;

  if(firsttime){
    
    firsttime=0;
    g=malloc(sizeof(double)*pos.nAtoms*3);
    h=malloc(sizeof(double)*pos.nAtoms*3);
    if(g==NULL|| h==NULL){
      printf("out of memory\n");
      exit(2);
    }
    initCgPar();	 
  }

  numOfCalls++;
  checkOptim=(numOfCalls%cgPar.part.optim.checkInterval)==0; /*this is slightly wierd. numOfCalls is saved over 
																					different minimizations => parts are increased at
																					different steps within minimzations*/
  
    
  fp=(*partpot)(par,pos,part,atomsInPart,1,1,0,cgPar.part.enlTreshold); /* the parts are now initialized.. */
  xi=pos.xa;
  potVal=fp;

  for (j=0;j<3*pos.nAtoms;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }

  
  for (its=0;its<ITMAX;its++) {
    partAtFinalSize=((*atomsInPart)==pos.nAtoms || its>cgPar.part.enlargements*cgPar.part.enlargeInterval);
    
    linminPart(&pos,h,par,&potVal,cgPar.part.allowAdaptive,part,(*atomsInPart),partpot);
    
    if(!restPotCalculated && partAtFinalSize){ /*this is activated first time that the part is at its final size*/
      restPotCalculated=1;
      restPot=(*partpot)(par,pos,part,atomsInPart,1,0,1,cgPar.part.enlTreshold)-fp; /*this is how much energy is in the bonds & angles not calculated now,no part of part */
      if(DEBUG) printf("partSize %d\n",(*atomsInPart));
      if(meanPartSize<0)
		  meanPartSize=(*atomsInPart);
		else
		  meanPartSize=0.99*meanPartSize+0.01*(*atomsInPart);
    }
    
    if (partAtFinalSize && 2.0*fabs(potVal-fp) <= cgPar.part.tolerance*(fabs(potVal)+fabs(fp)+EPS)) {
      break;	
    }
#ifdef PARALLEL
    { /*this block aborts this minimization if some oter process has an accepted step */
      int flag;
      MPI_Status status;
      MPI_Iprobe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&flag,&status); /*probe if somebody has sent a message, i.e. if a step has been accepted*/
      if(flag){
	parallelAbortion=1;
	break;
      }
    }
#endif
    
    
    fp=(*partpot)(par,pos,part,atomsInPart,1,((its%cgPar.part.enlargeInterval==0) && !partAtFinalSize ),0,cgPar.part.enlTreshold); 
    /*xi=pos.xa updated*/
     
    totForce=dgg=gg=0.0;
    for(jps=0;jps<(*atomsInPart);jps++) for(d=0;d<3;d++){
      j=part[jps]+d*pos.nAtoms;	
      totForce+=xi[j]*xi[j];
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
	 
	 
	 
    if (gg == 0.0) 
      break;
	 
    if(cgPar.part.optim.useOptim && partAtFinalSize  && (restPot+fp)-cgPar.part.optim.c*totForce>accPotVal){
      if(checkOptim && !checkOptimAborted){
		  checkOptimAborted=1;		  /* ok tried to break, used later on to see if c value was alright or not */
		  checkOptimPotential=restPot+fp;
		  checkOptimForce=totForce;
      }
      else{
		  if(DEBUG) printf("partCg aborted on step %d\n",its);
		  break;
      }
    }
	 
    gam=dgg/gg;
    for(jps=0;jps<(*atomsInPart);jps++) for(d=0;d<3;d++){
      j=part[jps]+d*pos.nAtoms;	
      g[j] = -xi[j];
      h[j]=g[j]+gam*h[j];
    }
    
  }
  
  /* 
     the final energy calculation done here have to be with fullenergy otherwise pos.enPot[] might be polluted with 
     false values near the part border.
  */   
  
  potVal=(*partpot)(par,pos,part,atomsInPart,1,0,1,cgPar.part.enlTreshold); 
  if(parallelAbortion)
    potVal=1.0e30; /*do not accept*/

  if(!parallelAbortion && cgPar.part.optim.useOptim && checkOptim && checkOptimAborted){
    double cGuess=1.1*(checkOptimPotential-potVal)/checkOptimForce; /*corrects c value to a more realistic one.. */
    if(cgPar.part.optim.c<cGuess)
      cgPar.part.optim.c=cGuess;
    
    if(DEBUG) printf("cgpar_part_optim_c: %g\n",cgPar.part.optim.c);
  }
  
  
  /*  printf("partcg:   part-full %g\n ",(potVal-pot(par,pos,0))/potVal);*/
	
  if(DEBUG) printf("%d steps in partCg, pot=%g\n",its,potVal/pos.nAtoms*TOEV);

  return potVal;
}  
 

/*This function first does a partial conj grad with optim, if it believes that this will probably be accepted at the step where the opim check is done it does a CG minimization of the whole system, otherwise it just returns a energy indicating that it is not accepted.*/
double initialPartCg(struct systemPos pos,struct parameters *par,int *part, int *atomsInPart,double accPotVal,	
		     double (*partpot)(struct parameters*,struct systemPos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold),
		     double (*pot)(struct parameters*,struct systemPos,int calcForces))
{
  int j,jps,d,its;
  static int firsttime=1;
  static int numOfCalls=0;
  
  double potVal;
  double gg,gam,fp,dgg;
  static double *g,*h,*xi;
  double totForce;
  int partAtFinalSize;
  
  int checkOptim=0;
  double checkOptimPotential;
  double checkOptimForce;
  int checkOptimAborted=0;

  if(firsttime){
    firsttime=0;
    g=malloc(sizeof(double)*pos.nAtoms*3);
    h=malloc(sizeof(double)*pos.nAtoms*3);
    initCgPar();
  }
  
  if(cgPar.initPart.doInitPart){
    numOfCalls++;
    checkOptim=(numOfCalls%cgPar.initPart.part.optim.checkInterval)==0;
	 
	 
    fp=(*partpot)(par,pos,part,atomsInPart,1,1,0,cgPar.initPart.part.enlTreshold); /* the parts are now initialized.. */
    xi=pos.xa;
    potVal=fp;
    for (j=0;j<3*pos.nAtoms;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j];
    }
	 
	 
    for (its=1;its<cgPar.initPart.part.enlargements*cgPar.initPart.part.enlargeInterval;its++) {
      partAtFinalSize=((*atomsInPart)==pos.nAtoms);

      linminPart(&pos,h,par,&potVal,cgPar.initPart.allowAdaptive,part,(*atomsInPart),partpot);
      fp=(*partpot)(par,pos,part,atomsInPart,1,((its%cgPar.initPart.part.enlargeInterval==0) &&!partAtFinalSize),0,cgPar.initPart.part.enlTreshold); 
      /*we already have the potential, we are interested in the new forces and the updating of the part*/
	 
      dgg=gg=0.0;
      for(jps=0;jps<(*atomsInPart);jps++) for(d=0;d<3;d++){
		  j=part[jps]+d*pos.nAtoms;	
		  gg += g[j]*g[j];
		  dgg += (xi[j]+g[j])*xi[j];
      }
      
      if (gg == 0.0) 
		  break;
		
      gam=dgg/gg;
	  
      for(jps=0;jps<(*atomsInPart);jps++) for(d=0;d<3;d++){
		  j=part[jps]+d*pos.nAtoms;	
		  g[j] = -xi[j];
		  h[j]=g[j]+gam*h[j];
      }
	 
    }
	 
    if(DEBUG) printf("initialpartSize  %d\n",(*atomsInPart));
	 if(meanInitPartSize<0) /*first time*/
		meanInitPartSize=(*atomsInPart);
	 else /*otherwise slowly changing*/
		meanInitPartSize=0.99*meanInitPartSize+0.01*(*atomsInPart);
	 
    potVal=(*partpot)(par,pos,part,atomsInPart,1,0,1,cgPar.initPart.part.enlTreshold); /*this is current total energy */
    /*printf("initPart: part-full %g\n ",(potVal-pot(par,pos,0))/potVal);*/ /*check that the total energy is correct*/
	 
    if(cgPar.initPart.part.optim.useOptim){
      totForce=0;
      
      for(jps=0;jps<(*atomsInPart);jps++) for(d=0;d<3;d++){
		  j=part[jps]+d*pos.nAtoms;	
		  totForce+=pos.xa[j]*pos.xa[j];
      }
		
      if(potVal-cgPar.initPart.part.optim.c*totForce>accPotVal){  /*if true then the optimization says that this step should be aborted, it will probably not be accepted anyway */
		  if(checkOptim){
			 checkOptimAborted=1;		  /* ok tried to break, used later on to see if c value was alright or not */
			 checkOptimPotential=potVal;
			 checkOptimForce=totForce;
		  }
		  else{
			 if(DEBUG) printf("initialPartCg aborted on step %d\n",its);
			 return potVal;
		  }
      }
    }
    
  }
  
#ifdef PARALLEL
  { /*this block aborts this minimization if some oter process has an accepted step */
    int flag=0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&flag,&status); /*probe if somebody has sent a message, i.e. if a step has been accepted*/
    if(flag){
      potVal=1e30; /* do not accept this step*/
      return potVal;
    }
  }
#endif
  
  if(cgPar.initPart.doPartCg)
    potVal=partCg(pos,par,part,atomsInPart,accPotVal,partpot,pot);
  if(cgPar.initPart.doFullCg)
    potVal=fullCg(pos,par,1,0,accPotVal,pot); /*minimize energy*/ /*we took away the isHarmonic part, perhaps its better not to have it.. */
  
  
  if( cgPar.initPart.part.optim.useOptim && checkOptim && checkOptimAborted){
    double cGuess=1.1*(checkOptimPotential-potVal)/checkOptimForce; /*corrects c value to a more realistic one.. */
    if(cgPar.initPart.part.optim.c<cGuess)
      cgPar.initPart.part.optim.c=cGuess;
    
    if(DEBUG) printf("cgpar_initialpart_optim_c: %g\n",cgPar.initPart.part.optim.c);
  }
  
  
  if(DEBUG) printf("%d steps in initialPartCg\n",its);
  
  return potVal;
}  


