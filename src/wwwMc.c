/*
BOMC - Monte Carlo based generator of defect-free amorphous samples of Si and SiO2

If you publish research done using Bomc then you are kindly asked to cite: 
 S. von Alfthan, A. Kuronen, and K. Kaski, Phys. Rev. B 68, 073203 (2003).


Copyright Sebastian von Alfthan  (galfthan at iki dot fi) 2002, 2003, 2014.


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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include "shared.h"
#include "fileio.h"
#include "potentialWww.h"
#include "mcSubs.h"
#include "miscSubs.h"
#include "conjgrad.h"
#include "parallel.h"
#include "randFunc.h"
#include "bondList.h"                     
#include "rings.h"  
#include "initialize.h"  
#include "wwwMc.h"

#define MAXERRORCODES 10

struct simuRes{
  int trialMoves[6];
  int accTrialMoves[6];
  int accMoves[6];
  int totAccTrialMoves;
  
};

static struct simuRes sres;       /*contains values for this processor, updated each trialstep */
static struct simuRes sresTotal; /*this is updated after each accepted step (with total values,summed over all processors)*/

void makeMove(struct systemPos pos,struct parameters *par,int accepted,double potVal,double *oldPotVal);

int prepareInitialSystem(struct systemPos pos,struct parameters *par,struct simuPar *simuPar){
  int rank;
  double potVal;
  int i;
  int errorCodes[MAXERRORCODES];
  int part[100]; /*small temporary part... */
  int atomsInPart;

  rank=getRank();
  if(rank==0){
    printf("Minimizing energy of original system from %g...",wwwPot(par,pos,0)*TOEV/pos.nAtoms);
    fflush(stdout);
  }
  
  potVal=fullCg(pos,par,0,0,1e10,wwwPot); /*optimization not allowed, accPotVal has a arbitrary value */
  enforcePeriodicity(pos,par);
  if(rank==0)
    printf("ready. System minimized to  pot %g\n",potVal*TOEV/pos.nAtoms);
  for(i=0;i<MAXERRORCODES;i++)
    errorCodes[i]=0;
  
  /*first make random bondswitches and accept them all, this is used to randomize the structre*/
  if(simuPar->randIter>0){
    if(rank==0) printf("Randomizing System\n");
    for(i=0;i<simuPar->randIter;i++)
      errorCodes[mcBondSwitch(&pos,part,&atomsInPart)]++;
    
    potVal=fullCg(pos,par,0,0,1e10,wwwPot); /*optimization not allowed, accPotVal has a arbitrary value */
    enforcePeriodicity(pos,par);
    
    if(rank==0){
      printf("Randomizing part ready.  \n");
      for(i=0;i<MAXERRORCODES;i++)
	if(errorCodes[i]>0)
	  printf("switch err: %d %d\n",i,errorCodes[i]);
      printf("After %d randsteps, pot:%geV\n",simuPar->randIter,potVal*TOEV/pos.nAtoms);  
    }


  }
  backupState(pos);
  return 0;
}


int mcSimulation(struct systemPos pos,struct parameters *par,struct simuPar *simuPar,struct fileIo *fio){
  int *part;
  int atomsInPart;
  int rank;
  double rNum;
  double potVal,oldPotVal;
  double oldVol;
  int i;
  time_t startTime;		
  char buffer[512];
  static int prevWriteStep=0; /*at what step did we last time write outthe whole system*/
  
  rank=getRank();
  /* initialize these arrays */
  for(i=0;i<6;i++){ 
    sres.trialMoves[i]=0;				  /* how many trial moves */
    sres.accTrialMoves[i]=0;			  /* hom many of the trial moves we are abel to try */
    sres.accMoves[i]=0;				  /* how many acepted moves */
  }
  sres.totAccTrialMoves=0;

  part=malloc(pos.nAtoms*sizeof(int));

  /*update pots*/
  potVal=wwwPot(par,pos,0);
  /*now the potential and position of the atoms have to be saved so that makeMove can use the correct backups*/
  oldPotVal=potVal;
  backupState(pos);

  reinitWwwPotAfterStep(par,pos,1); /*make sure that wwwPot properly initialized*/

  startTime=time(NULL);
  while((par->kT=getKT(sres.totAccTrialMoves,pos.nAtoms))>=0.0){ /* while we are in the annealing schedule */
    int error;
    int mcStep;
    /* let chose mcStep */
    if(par->kT>simuPar->volStartkT || randNum(0)>simuPar->volProb){ /* bond switch or vomulme step */
      rNum=randNum(0);
      if(rNum<simuPar->bondSwitchProb)
	mcStep=0;				  /* switch */
      else if(rNum<simuPar->bondSwitchProb+simuPar->bondBreakProb)
	mcStep=1;				  /* break */
      else if(rNum<simuPar->bondSwitchProb+simuPar->bondBreakProb+simuPar->bondCreateProb)
	mcStep=2;				  /* create */
      else if(rNum<simuPar->bondSwitchProb+simuPar->bondBreakProb+simuPar->bondCreateProb+simuPar->bondDiffuseProb)
	mcStep=3;					  /* diffuse */
      else
	mcStep=5; /*oxygen diffusion*/
    }
    else
      mcStep=4;					  /* volume */
	 
    
    switch(mcStep){ 								  /* step chosen, make it */
    case 0:
      error=mcBondSwitch(&pos,part,&atomsInPart);
      break;
    case 1:
      error=mcBondBreak(&pos,part,&atomsInPart);
      break;
    case 2:	
      error=mcBondCreate(&pos,par,part,&atomsInPart);
      break;
    case 3:
      error=mcBondDiffuse(&pos,part,&atomsInPart);
      break;
    case 4:
      oldVol=par->volume;
      mcVolume(pos,par,0);
      error=0;						  /*volume  always succeeds. */
      break;
    case 5:
      error=mcOxygenDiffuse(&pos,part,&atomsInPart);
      break;
    }
    
    reinitWwwAfterBondChange(par,pos);  /*reinits potential after bond change */
    sres.trialMoves[mcStep]++;
    
    if(error==0){				  /* if we made an accepted trial move */
      double accPotVal;
      sres.totAccTrialMoves++;
      sres.accTrialMoves[mcStep]++;
      
      /* minimize energy of new state, depends on if it was a vilume step or not */
      if(mcStep==4){				  /* volume step */
		  accPotVal=oldPotVal-par->pressure*(par->volume-oldVol)-par->kT*log(randNum(0))+par->kT*(pos.nAtoms+1)*log(par->volume/oldVol) ;
		  potVal=fullCg(pos,par,1,0,accPotVal,wwwPot);
      }
      else{
	accPotVal=oldPotVal-par->kT*log(randNum(0));
	potVal=initialPartCg(pos,par,part,&atomsInPart,accPotVal,partwwwPot,wwwPot); 
      }
      
      if(potVal<accPotVal){ /*move accepted*/
	sres.accMoves[mcStep]++;
	makeMove(pos,par,1,potVal,&oldPotVal); 
      }
      
      else{							  /* move not accepted */
	if(mcStep==4) mcVolume(pos,par,1); /*restore old Vol*/
	makeMove(pos,par,0,potVal,&oldPotVal); /*not accepted*/
	potVal=oldPotVal;
      }
		
      /* Write out some results and write out the system */
      if(sresTotal.totAccTrialMoves>simuPar->resultPeriod+prevWriteStep && getRank()==0){ 
	double extime;
	int days,hours,minutes,seconds;
	prevWriteStep=sresTotal.totAccTrialMoves;
	sprintf(buffer,"steps %d kT %g pot %g density %g",sresTotal.totAccTrialMoves, par->kT*TOEV,potVal*TOEV/pos.nAtoms,pos.nAtoms/par->volume);
	extime=difftime(time(NULL),startTime)+1; /* +1 to avoid problems if one gets here in under one seconds, avoid a divide by zero */
	convertSeconds(extime,&days,&hours,&minutes,&seconds);
	sprintf(buffer,"%s elapsed time %dd %dh %dm",buffer,days,hours,minutes);
	enforcePeriodicity(pos,par);
	writeSystem(fio,&pos,par,buffer);
	getCgStatus(buffer);
	if(getRank()==0) printf("cgStatus: %s\n",buffer);
	if(rank==0) fflush(stdout);			  /* annoying when it doesn't write to logfile.. */
      }
    }
    
  }
  return 0;
}



/* two verions of makeMove, parallel and serial*/

#ifndef PARALLEL
/*this is the serial version*/
void makeMove(struct systemPos pos,struct parameters *par,int accepted,double potVal,double *oldPotVal){
  int i;
  static int accSteps=0;
  static int steps=0;
  
  int dBonds=0;
  struct bondList *blist;
  getBondList(&blist);
  
  for(i=0;i<6;i++){
    sresTotal.trialMoves[i]=sres.trialMoves[i];
    sresTotal.accTrialMoves[i]=sres.accTrialMoves[i];
    sresTotal.accMoves[i]=sres.accMoves[i];
  }
  sresTotal.totAccTrialMoves=sres.totAccTrialMoves;


  steps++;  
  if(accepted){
    for(i=0;i<pos.nAtoms;i++)
      dBonds+=blist->danglBonds[i];
    *oldPotVal=potVal;
    accSteps++;
    
    enforcePeriodicity(pos,par);
    backupState(pos);
    
    printf("%6d accsteps  tot: %d/%d= %.2g  s %d/%d= %.2g bb %d/%d= %.2g bc %d/%d= %.2g bd %d/%d= %.2g v %d/%d= %.2g od %d/%d= %.2g\n" 
	   ,accSteps,accSteps,steps,(double)accSteps/steps
	   ,sresTotal.accMoves[0],sresTotal.accTrialMoves[0],sresTotal.accMoves[0]/(1e-10+sresTotal.accTrialMoves[0])
	   ,sresTotal.accMoves[1],sresTotal.accTrialMoves[1],sresTotal.accMoves[1]/(1e-10+sresTotal.accTrialMoves[1])
	   ,sresTotal.accMoves[2],sresTotal.accTrialMoves[2],sresTotal.accMoves[2]/(1e-10+sresTotal.accTrialMoves[2])
	   ,sresTotal.accMoves[3],sresTotal.accTrialMoves[3],sresTotal.accMoves[3]/(1e-10+sresTotal.accTrialMoves[3])
	   ,sresTotal.accMoves[4],sresTotal.accTrialMoves[4],sresTotal.accMoves[4]/(1e-10+sresTotal.accTrialMoves[4])
	   ,sresTotal.accMoves[5],sresTotal.accTrialMoves[5],sresTotal.accMoves[5]/(1e-10+sresTotal.accTrialMoves[5]));
    
    printf("%6d dens: %g pot: %.10g kT: %gev %g db  \n"
	   ,accSteps,(double)pos.nAtoms/par->volume,potVal*TOEV/pos.nAtoms
	   ,par->kT*TOEV,(double)dBonds/pos.nAtoms ); 	
    fflush(stdout);
    reinitWwwAfterBondChange(par,pos);
  }
  else /*no accepted move..*/
    restoreState(pos);

  reinitWwwPotAfterStep(par,pos,accepted);
}
#endif

#ifdef PARALLEL
/*this is parallel version*/
void makeMove(struct systemPos pos,struct parameters *par,int accepted,double potVal,double *oldPotVal){
  MPI_Status status;
  MPI_Status status2;
  int flag;
  int rank,totp;
  static int *accVect;
  static int *accVectRes;
  int chosenRank;
  int i,j;
  static int firsttime=1;
  static int accSteps=0;
  static int steps=0;
  int totSteps;
  static double waitTime=0.0;
  static double totTime=0.0;
  struct bondList *blist;
  getBondList(&blist);
  steps++;
  rank=getRank();
  totp=getCommSize();
  
  if(firsttime){
    accVect=malloc(sizeof(int)*totp); /*this vector will tell us which processors have made an accepted step*/
    accVectRes=malloc(sizeof(int)*totp); /*this is the reuslt when we sum all the accVEcts in all processors*/
    firsttime=0;
  }
   

  MPI_Iprobe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&flag,&status); /*probe if somebody has sent a message, i.e. if a step has been accepted*/
  
  if(flag){
    int dummy;
    MPI_Recv(&dummy,1,MPI_INT,status.MPI_SOURCE,2,MPI_COMM_WORLD,&status2); /*clear away recvs*/
  }
  else if(accepted){ /*accepted and no one else has yet been acepted, make a send so that they also stop*/
    int dummy=1,i;
    for(i=0;i<totp;i++){
      if(i==rank) continue;
      MPI_Send(&dummy,1,MPI_INT,i,2,MPI_COMM_WORLD);
    }
  }
  
  if(flag || accepted) { /*somebody has been acceepted*/
    for(i=0;i<totp;i++){
      accVect[i]=0;
      accVectRes[i]=0;
    }
    accVect[rank]=accepted;
    
    MPI_Reduce(accVect,accVectRes,totp,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); /*sum the total acceptance vector at rank=0*/
    

    if(rank==0){ /*rank=0 chooses the accepted step*/
      j=0;      
      /*first we change the format of accVect if totp=4 and 2 and 3 had accepted steps then
	accvect = 0 1 1 0  => 2 3 - - */
      for(i=0;i<totp;i++)
	if(accVectRes[i]){
	  accVectRes[j]=i;
	  j++;
	}
      if(j!=0)
	chosenRank=accVectRes[(int)(randNum(2)*j)]; /*we take on random one of the chosen events*/
      else
	chosenRank=-1;
    }
    
    MPI_Bcast(&chosenRank,1,MPI_INT,0,MPI_COMM_WORLD); /*send to everybody who was chosen*/

    if(chosenRank==-1 ){ /*some error happened (bugs...) when we chose the processor, it actuallly wasnt accepted*/
      printf("error: None accepted rank:%d chosenrank=%d\n",rank,chosenRank);
      restoreState(pos);
      reinitWwwPotAfterStep(par,pos,0); 
    }
    
    else{
      accSteps++;
      MPI_Bcast(par->box,3,MPI_DOUBLE,chosenRank,MPI_COMM_WORLD);
      MPI_Bcast(pos.x,3*pos.nAtoms,MPI_DOUBLE,chosenRank,MPI_COMM_WORLD);
      MPI_Bcast(pos.enPot,pos.nAtoms,MPI_DOUBLE,chosenRank,MPI_COMM_WORLD);
      MPI_Bcast(blist->head,blist->arraySize,MPI_INT,chosenRank,MPI_COMM_WORLD);
      MPI_Bcast(&potVal,1,MPI_DOUBLE,chosenRank,MPI_COMM_WORLD);
      *oldPotVal=potVal;
      
      setSystemSize(par);
      enforcePeriodicity(pos,par);
      backupState(pos);
      reinitWwwAfterBondChange(par,pos);
      reinitWwwPotAfterStep(par,pos,1);
      
      MPI_Reduce(&steps,&totSteps,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 
      MPI_Reduce(&sres.totAccTrialMoves,&sresTotal.totAccTrialMoves,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 
      MPI_Reduce(&(sres.trialMoves[0]),&(sresTotal.trialMoves[0]),6,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 
      MPI_Reduce(&(sres.accTrialMoves[0]),&(sresTotal.accTrialMoves[0]),6,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 
      MPI_Reduce(&(sres.accMoves[0]),&(sresTotal.accMoves[0]),6,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 

      
      if(rank==0) {
	int dBonds=0;
	for(i=0;i<pos.nAtoms;i++)
	  dBonds+=blist->danglBonds[i];
	printf("%6d accsteps  tot: %d/%d= %.2g  s %d/%d= %.2g bb %d/%d= %.2g bc %d/%d= %.2g bd %d/%d= %.2g v %d/%d= %.2g od %d/%d= %.2g\n" 
	       ,accSteps,accSteps,totSteps,(double)accSteps/totSteps
	       ,sresTotal.accMoves[0],sresTotal.accTrialMoves[0],sresTotal.accMoves[0]/(1e-10+sresTotal.accTrialMoves[0])
	       ,sresTotal.accMoves[1],sresTotal.accTrialMoves[1],sresTotal.accMoves[1]/(1e-10+sresTotal.accTrialMoves[1])
	       ,sresTotal.accMoves[2],sresTotal.accTrialMoves[2],sresTotal.accMoves[2]/(1e-10+sresTotal.accTrialMoves[2])
	       ,sresTotal.accMoves[3],sresTotal.accTrialMoves[3],sresTotal.accMoves[3]/(1e-10+sresTotal.accTrialMoves[3])
	       ,sresTotal.accMoves[4],sresTotal.accTrialMoves[4],sresTotal.accMoves[4]/(1e-10+sresTotal.accTrialMoves[4])
	       ,sresTotal.accMoves[5],sresTotal.accTrialMoves[5],sresTotal.accMoves[5]/(1e-10+sresTotal.accTrialMoves[5]));
    	
	printf("%6d dens: %g pot: %.10g kT: %gev %g db r: %d acc: %d/%d\n"
	       ,accSteps,(double)pos.nAtoms/par->volume,potVal*TOEV/pos.nAtoms
	       ,par->kT*TOEV,(double)dBonds/pos.nAtoms,chosenRank,j,totp); 	
	fflush(stdout);
      }
    }
  }
  
  /*no accepted move..*/
  else {
    restoreState(pos);
    reinitWwwPotAfterStep(par,pos,0);
  }
  
}
#endif
