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



/*! 
  \file  potentialWww.c
  \brief Module caclulating potential energy and forces of system
  
  This module calculates the potential energy of the system. There are two 5 different part of the potential.
  -Normal Keating
  -Simplified Keating
  -Repulsive 
  -Dangling bonds
  -Suboxide penalty


  The Normal Keating potential is as follows:
  
  V=sum_{bonds} V_{ij} +sum_{angles} V_{ijk}

  where
  
  V_{ij} =K_s*[r0_{ij}^2-r_{ij}^2]^2;
  V_{ijk}=K_b*[r_i \ldot r_j - cos(\theta_0)*r0_{ij}*r0_{ik}]^2;
  

  The Simpilfied Keating potential is as follows:
  
  V=sum_{bonds} V_{ij} +sum_{angles} V_{ijk}

  where
  
  V_{ij} =0.5*K_s*[r0_{ij}-r_{ij}]^2;
  V_{ijk}=0.5*K_b*[cos(\theta_{ijk}) - cos(\theta_0)]^2;
  
 

*/
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include "shared.h"
#include "fileio.h"
#include "ngbrs.h"
#include "potentialWww.h"
#include "bondList.h"
#include "parallel.h"
#include "rings.h"
#include "tabulatedValues.h"
static int siType=1,oType=0;

struct wwwPotParameters{
  double repR0;
  double repR02;
  double repK;
  double rcut,rcut2;
  
  double r0[2][2];
  double r02[2][2];
  double cos0[2];
  double kStretch[2][2];
  double kBend[2][2][2];
  
  double subOxidePenalty[5];
  
  double ringPenalty3;
  double ringPenalty2;
  
  
  struct valueTable ooPottable;
  struct valueTable soPottable;
  struct valueTable ssPottable;
  int includeRepPot;
  int includeKeatingPot;
  int includeSuboxidePot;
  int includeDanglbondPot;
  int includeRingPot;
  int includeSnnRep;
  double ringPenaltyEnergy;

  double stix_D;
  double stix_beta;
  double stix_r0;
  double stix_galpha;
  double stix_alpha0;
  double stix_gL;
  double stix_L0;
  double stix_A;
  double stix_b;
  double stix_rc;
  double stix_gamma;
  struct valueTable stix_repTable;
};


#define KBEND(i,j,k) (wwwPar.kBend[pos.aType[i]][pos.aType[j]][pos.aType[k]])
#define KSTRETCH(i,j) (wwwPar.kStretch[pos.aType[i]][pos.aType[j]])
#define COS0(i) (wwwPar.cos0[pos.aType[i]])
#define R0(i,j) (wwwPar.r0[pos.aType[i]][pos.aType[j]])
#define R02(i,j) (wwwPar.r02[pos.aType[i]][pos.aType[j]])


static struct wwwPotParameters wwwPar;
static struct ngbrData nd;
static struct ngbrData ndLastStep;
static int *prevInteractionHead;
static int *prevInteractionList;
static double prevRingPenaltyEnergy;

/* prototypes of internal functions */

static double partLoop(struct parameters *par,struct systemPos pos,int *partSkin,int atomsInPartSkin,int calcForces,
							  double (*innerLoop)(int i,struct parameters *par,struct systemPos pos,int calcForces));
static double allLoop(struct parameters *par,struct systemPos pos,int calcForces,
							 double (*innerLoop)(int i,struct parameters *par,struct systemPos pos,int calcForces));

static double keatingInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces);
static double simpKeatingInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces);
static double repulsiveInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces);
static double suboxideInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces);
static double danglBondsInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces);
static double tabulatedRepulsiveInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces);

static void updateParts(struct parameters *par,struct systemPos pos,int *part,int *partSkin, int *atomsInPart,int *atomsInPartSkin,int includePrevInteractionsInSkin,int updatePart);

static void putInteractingInTable(int atom,struct parameters *par,struct systemPos *pos,struct bondList *blist,int *table);
static void putPrevInteractingInTable(int atom,struct parameters *par,struct systemPos *pos,int *table);
static void putInteractingInPrevList(int atom,struct parameters *par,struct systemPos *pos,struct bondList *blist,  int *posInPrevList);
double ringPenaltyEnergy(struct parameters *par,struct systemPos pos);
static double stixInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces);


/*include the innerloop funxctions and the other functions that calculate the real potentials */
#include "forceFunctions.c"


void initwwwPot(struct parameters *par,struct systemPos pos){
  int error,i,j,k;
  char filename[20];
  struct bondList *blist;
  getBondList(&blist);
  
  strcpy(filename,"wwwPot.par");
  openParameterFile(filename);
  error=0;
  
  error+=getParValue("includeKeatingPot",&wwwPar.includeKeatingPot,"%d");
  error+=getParValue("includeRepPot",&wwwPar.includeRepPot,"%d");
  error+=getParValue("includeSnnRep",&wwwPar.includeSnnRep,"%d");
  error+=getParValue("includeSuboxidePot",&wwwPar.includeSuboxidePot,"%d");
  error+=getParValue("includeDanglbondPot",&wwwPar.includeDanglbondPot,"%d");
  error+=getParValue("includeRingPot",&wwwPar.includeRingPot,"%d");
  
  if(wwwPar.includeKeatingPot==1){/*original keatin */
    error+=getParValue("okp1_kStretch_SiSi",&wwwPar.kStretch[1][1],"%lf");
    error+=getParValue("okp1_kStretch_SiO",&wwwPar.kStretch[1][0],"%lf");
    error+=getParValue("okp1_kBend_SiSiSi",&wwwPar.kBend[1][1][1],"%lf");
    error+=getParValue("okp1_kBend_SiOSi",&wwwPar.kBend[1][0][1],"%lf");
    error+=getParValue("okp1_kBend_OSiO",&wwwPar.kBend[0][1][0],"%lf");
    error+=getParValue("okp1_r0_SiSi",&wwwPar.r0[1][1],"%lf");
    error+=getParValue("okp1_r0_SiO",&wwwPar.r0[1][0],"%lf");
    error+=getParValue("okp1_cos0_Si",&wwwPar.cos0[1],"%lf");
    error+=getParValue("okp1_cos0_O",&wwwPar.cos0[0],"%lf");
  }
  else { /*simplified keatin */
    error+=getParValue("skp2_kStretch_SiSi",&wwwPar.kStretch[1][1],"%lf");
    error+=getParValue("skp2_kStretch_SiO",&wwwPar.kStretch[1][0],"%lf");
    error+=getParValue("skp2_kBend_SiSiSi",&wwwPar.kBend[1][1][1],"%lf");
    error+=getParValue("skp2_kBend_SiOSi",&wwwPar.kBend[1][0][1],"%lf");
    error+=getParValue("skp2_kBend_OSiO",&wwwPar.kBend[0][1][0],"%lf");
    error+=getParValue("skp2_r0_SiSi",&wwwPar.r0[1][1],"%lf");
    error+=getParValue("skp2_r0_SiO",&wwwPar.r0[1][0],"%lf");
    error+=getParValue("skp2_cos0_Si",&wwwPar.cos0[1],"%lf");
    error+=getParValue("skp2_cos0_O",&wwwPar.cos0[0],"%lf");
	 
  }
  if(wwwPar.includeKeatingPot==3){/*stixrudes potential, no keating at all...*/
	 error+=getParValue("stix_D",&wwwPar.stix_D,"%lf");
	 error+=getParValue("stix_beta",&wwwPar.stix_beta,"%lf");
	 error+=getParValue("stix_r0",&wwwPar.stix_r0,"%lf");
	 error+=getParValue("stix_galpha",&wwwPar.stix_galpha,"%lf");
	 error+=getParValue("stix_alpha0",&wwwPar.stix_alpha0,"%lf");
	 error+=getParValue("stix_gL",&wwwPar.stix_gL,"%lf");
	 error+=getParValue("stix_L0",&wwwPar.stix_L0,"%lf");
	 error+=getParValue("stix_A",&wwwPar.stix_A,"%lf");
	 error+=getParValue("stix_b",&wwwPar.stix_b,"%lf");
	 error+=getParValue("stix_rc",&wwwPar.stix_rc,"%lf");
	 error+=getParValue("stix_gamma",&wwwPar.stix_gamma,"%lf");
 
  }

  
  
  error+=getParValue("repulsive_r0",&wwwPar.repR0,"%lf");
  error+=getParValue("repulsive_k",&wwwPar.repK,"%lf");
  
  error+=getParValue("ringPenalty3",&wwwPar.ringPenalty3,"%lf");
  error+=getParValue("ringPenalty2",&wwwPar.ringPenalty2,"%lf");

  error+=getParValue("SubOxidePenalty+0",&wwwPar.subOxidePenalty[0],"%lf");
  error+=getParValue("SubOxidePenalty+1",&wwwPar.subOxidePenalty[1],"%lf");
  error+=getParValue("SubOxidePenalty+2",&wwwPar.subOxidePenalty[2],"%lf");
  error+=getParValue("SubOxidePenalty+3",&wwwPar.subOxidePenalty[3],"%lf");
  error+=getParValue("SubOxidePenalty+4",&wwwPar.subOxidePenalty[4],"%lf");
  closeParameterFile();
  
  for(i=0;i<2;i++) 
    for(j=0;j<2;j++){
      wwwPar.kStretch[i][j]*=0.5; 
      for(k=0;k<2;k++)
		  wwwPar.kBend[i][j][k]*=0.5;
    }  
  wwwPar.repK*=0.5;
  
  wwwPar.kBend[1][1][0]=sqrt(wwwPar.kBend[1][1][1]*wwwPar.kBend[0][1][0]);
  wwwPar.kBend[0][1][1]=wwwPar.kBend[1][1][0];
  
  
  /*make the matrices symmetric*/
  wwwPar.kStretch[0][1]=wwwPar.kStretch[1][0];
  wwwPar.r0[0][1]=wwwPar.r0[1][0];

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      wwwPar.r02[i][j]=POW2(wwwPar.r0[i][j]);
  
  /*change the units*/
  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      wwwPar.kStretch[i][j]*=EV;
      for(k=0;k<2;k++)
		  wwwPar.kBend[i][j][k]*=EV;
    }
  for(i=0;i<5;i++)
    wwwPar.subOxidePenalty[i]*=EV;

  wwwPar.ringPenalty3*=EV;  
  wwwPar.ringPenalty2*=EV;  
  wwwPar.repK*=EV;
  wwwPar.stix_D*=EV;
  wwwPar.stix_A*=EV;
  wwwPar.stix_galpha*=EV;
  wwwPar.stix_gL*=EV;
  if(wwwPar.includeKeatingPot==3)
	 readTable("stixrep.dat",&wwwPar.stix_repTable);

  wwwPar.rcut=wwwPar.repR0;
  wwwPar.rcut2=POW2(wwwPar.rcut);
  wwwPar.repR02=POW2(wwwPar.repR0);
  
  if(wwwPar.includeSnnRep){
	 readTable("ooPottable.dat",&wwwPar.ooPottable);
	 readTable("ssPottable.dat",&wwwPar.ssPottable);
	 readTable("soPottable.dat",&wwwPar.soPottable);
  }
  if(error!=0){
    if(getRank()==0) printf("could not read the required parameters from file %s\n",filename);
    exit(0);
  }

  if(wwwPar.includeRepPot && wwwPar.includeKeatingPot==3){
	 printf("stix and reppot cannot be on at the same time");
	 exit(0);
  }

  /*init ngbrs */
  initNgbrs(&nd,pos,par,wwwPar.rcut);
  initNgbrs(&ndLastStep,pos,par,wwwPar.rcut);
  
  updateNgbrs(&nd,pos,par);
  if(wwwPar.includeRepPot)
	 nonNnNgbrs(&nd,pos.nAtoms); 
  
  
  prevInteractionHead=malloc(sizeof(int)*pos.nAtoms);
  prevInteractionList=malloc(sizeof(int)*pos.nAtoms*(blist->nMax+nd.nMax)); /* if these lists are later on larger there might  occur a sigsegv */
  reinitWwwAfterBondChange(par, pos);
  reinitWwwPotAfterStep(par,pos,1); /* this updates some stuff.. (prevInteractionHead,prevInteractionList,ndLastStep) */
}


double wwwPot(struct parameters *par,struct systemPos pos,int calcForces){
  int i;
  double totPot=0;
  

  if(calcForces) /*this is not good if we want to have more than one pot at the same time!!!*/
    memset(pos.xa,0,sizeof(double)*3*pos.nAtoms);
  
  memset(pos.enPot,0,sizeof(double)*pos.nAtoms);
  
  /*calculate forces*/
  
  if(wwwPar.includeKeatingPot==1)  
    totPot+=allLoop(par,pos,calcForces,keatingInnerLoop);
  else if(wwwPar.includeKeatingPot==2)  /*simplified keatin */
    totPot+=allLoop(par,pos,calcForces,simpKeatingInnerLoop);
  else if(wwwPar.includeKeatingPot==3){  /*stix */
	 if(checkNgbrUpdateAll(&nd,pos,par))
      updateNgbrs(&nd,pos,par);
	 totPot+=allLoop(par,pos,calcForces,stixInnerLoop);
  }

  if(wwwPar.includeRepPot){
    if(checkNgbrUpdateAll(&nd,pos,par)){
      updateNgbrs(&nd,pos,par);
      nonNnNgbrs(&nd,pos.nAtoms); /*tkaes away first and scond ngbrs that are bonded*/
    }
    totPot+=allLoop(par,pos,calcForces,repulsiveInnerLoop);
  }

  if(wwwPar.includeSnnRep){
    totPot+=allLoop(par,pos,calcForces,tabulatedRepulsiveInnerLoop);
  }


  if(wwwPar.includeSuboxidePot)
    totPot+=allLoop(par,pos,calcForces,suboxideInnerLoop);
  if(wwwPar.includeDanglbondPot)
    totPot+=allLoop(par,pos,calcForces,danglBondsInnerLoop);
  if(wwwPar.includeRingPot)
    totPot+=wwwPar.ringPenaltyEnergy;

  if(calcForces)
    for(i=0;i<pos.nAtoms;i++)
      if(pos.aFlag[i]==FIXED_FLAG){
		  pos.xa[i]=0;
		  pos.ya[i]=0;
		  pos.za[i]=0;
      }
  
  return totPot;
}

/*return true (1) if atom is in part*/
int isInPart(int atom, int *part,int partSize){
  int k;
  for(k=0;k<partSize;k++)
    if(part[k]==atom)
      return 1;
  return 0;
}

double partwwwPot(struct parameters *par,struct systemPos pos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold){
  static int firsttime=1;
  static int *partSkin; /*part skin is larger than part, it also contains a skin so that the correct forces &pots can be calculated*/
  static int *partSkinTemp;
  int *partSkinBackup;
  static int atomsInPartSkin;
  static int atomsInPartSkinTemp;
  int atomsInPartSkinBackup;
  int i,j,ip;
  double partSkinPot=0.0;
  struct bondList *blist;
  getBondList(&blist);
  
  if(firsttime){
    partSkin=malloc(sizeof(int)*pos.nAtoms);
    partSkinTemp=malloc(sizeof(int)*pos.nAtoms);
    atomsInPartSkin=0;
    firsttime=0;
  }
  

  if(updatePart && returnRealPot)
    updateParts(par,pos,part,partSkin,atomsInPart,(&atomsInPartSkin),1,1); 
  else if(updatePart)
    updateParts(par,pos,part,partSkin,atomsInPart,(&atomsInPartSkin),0,1);
  /*  update partSkin so that it contains all atoms interacting, or that used  to interact with part, part unchanged, look at treshold... */
  else if(returnRealPot)
    updateParts(par,pos,part,partSkin,atomsInPart,(&atomsInPartSkin),1,0);


  if(returnRealPot){
    /*if part large it is actually faster just to calculate the whole energy of the system, this particular value is just a guess.. */
    if((double)*atomsInPart/pos.nAtoms>0.5) 
      return wwwPot(par,pos,calcForces);
    /* 
       after this 
       update partSkinTemp so that it contains all atoms interacting with partSkin
       partSkin uncchanged,
    */
	 
    updateParts(par,pos,partSkin,partSkinTemp,(&atomsInPartSkin),(&atomsInPartSkinTemp),0,0);
	 
    /* now we save the old partskin and use the larger partskin, this makes sure that the results are correct.. (the pair forces are then correct in old partskin) */
    partSkinBackup=partSkin;	  /* save pointer */
    atomsInPartSkinBackup=	atomsInPartSkin;  /* save size*/ 
    partSkin=partSkinTemp;		  /* put partskin to point to temporary list */
    atomsInPartSkin=atomsInPartSkinTemp;
  }

  /* One would think that it would be enough to put the forces at part to zero but this leads to some small errors... why? have no idea, there might be some bug somewhere? */

  if(calcForces)
    memset(pos.xa,0,sizeof(double)*3*pos.nAtoms);
  
  for(ip=0;ip<atomsInPartSkin;ip++){
    i=partSkin[ip];
    pos.enPot[i]=0.0;
  }
  
  
  /*calculate forces*/
  
  if(wwwPar.includeKeatingPot==1)  /*dont know why one would want to exclude this part but...*/
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,keatingInnerLoop);
 
  else if(wwwPar.includeKeatingPot==2 )  /*simplified keating */
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,simpKeatingInnerLoop);
  
  else if(wwwPar.includeKeatingPot==3) {
	 if(checkNgbrUpdatePart(&nd,pos,par,partSkin,atomsInPartSkin)){
		updateNgbrsPart(&nd,pos,par,partSkin,atomsInPartSkin);
	 }
	 partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,stixInnerLoop); 
  }
  
  if(wwwPar.includeRepPot){
    if(checkNgbrUpdatePart(&nd,pos,par,partSkin,atomsInPartSkin)){
      updateNgbrsPart(&nd,pos,par,partSkin,atomsInPartSkin);
      nonNnNgbrsPart(&nd,pos.nAtoms,partSkin,atomsInPartSkin); /*tkaes away first and scond ngbrs that are bonded*/
    }
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,repulsiveInnerLoop);
  }

  if(wwwPar.includeSnnRep){
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,tabulatedRepulsiveInnerLoop);
  }


  if(wwwPar.includeSuboxidePot)
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,suboxideInnerLoop);
  if(wwwPar.includeDanglbondPot)
    partSkinPot+=partLoop(par,pos,partSkin,atomsInPartSkin,calcForces,danglBondsInnerLoop);
  
  /*put forces to zero at atoms which are fixed*/
  for(ip=0;ip<atomsInPartSkin;ip++){
    i=partSkin[ip];
	 if(pos.aFlag[i]==FIXED_FLAG){
		pos.xa[i]=0;
		pos.ya[i]=0;
		pos.za[i]=0;
	 }
  }
  
  /*now we have to fix the pot values so that tey are correct for each atom in the whole system*/
  if(returnRealPot){
    double totPot=0;


    if(wwwPar.includeRingPot)
      totPot+=wwwPar.ringPenaltyEnergy; /*increase pot value, this is not included in enPot[] */ 

    partSkin=partSkinBackup;	  /* restore pointer */
    atomsInPartSkin=atomsInPartSkinBackup;  /* restore size*/ 
	 
    ip=0;	 
    for(i=0;i<pos.nAtoms;i++){
      j=partSkin[ip];
      if(i==j){						  /* i is in partskin */
		  if(ip<atomsInPartSkin)
			 ip++;						  /* next time we check against the following one in the list */
      }

      else
		  pos.enPot[i]=pos.prevEnPot[i]; /* outside partskin, near partskin the potentials are polluted, copy old correct vaalues here */
      totPot+=pos.enPot[i];
    }
    return totPot;
  }
  return partSkinPot;
}

static double partLoop(struct parameters *par,struct systemPos pos,int *partSkin,int atomsInPartSkin,int calcForces,
							  double (*innerLoop)(int i,struct parameters *par,struct systemPos pos,int calcForces)){
  double totPot=0;
  int ips;
  
  for(ips=0;ips<atomsInPartSkin;ips++)
    totPot+=(*innerLoop)(partSkin[ips],par,pos,calcForces);
  
  return totPot;
}

static double allLoop(struct parameters *par,struct systemPos pos,int calcForces,
							 double (*innerLoop)(int i,struct parameters *par,struct systemPos pos,int calcForces)){
  double totPot=0;
  int i;
  
  for(i=0;i<pos.nAtoms;i++)
    totPot+=(*innerLoop)(i,par,pos,calcForces);
  return totPot;
}


/*
  Description:
  This function puts all the atoms interacting with atom in table by using the same bondlist and ngbrList as the real potential functions in the module. Table is an integer array where table[j]=1 if it interacts with atom

*/

static void putInteractingInTable(int atom,struct parameters *par,struct systemPos *pos,struct bondList *blist,int *table){
  int j,bondListPosj,ngbrListPosj;
  double drx_ij,dry_ij,drz_ij;
  double rij2;
  double r02;
  
  /* r02 is the distance at which ngbr atoms are included, shoud this include skin ? probably not! */
  /*  r02=POW2(sqrt(wwwPar.repR02)+par->skin); */
  r02=wwwPar.repR02;
  
  /* first put the atom itself in the table */
  table[atom]=1;
  
  /* then the atoms it is bonded to */
  bondListPosj=blist->head[atom]; /*where the bonds start*/
  while(blist->list[bondListPosj]!=-1) /*get bonds of i, stops when the bonds od atom i+1 comes*/
    table[blist->list[bondListPosj++]]=1;
	 
	 
  /* and last the atoms it is ngbr with */
  if(wwwPar.includeRepPot){
    ngbrListPosj=nd.head[atom]; /*where the ngbrs start*/
    while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=nd.list[ngbrListPosj++]; /*j is the ngbr*/
      if(!table[j]){
		  drx_ij=pos->x[j]-pos->x[atom]; 
		  dry_ij=pos->y[j]-pos->y[atom]; 
		  drz_ij=pos->z[j]-pos->z[atom];	
		  PERIODIC(drx_ij,dry_ij,drz_ij);
		  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		  if( rij2 <r02 )
			 table[j]=1;
      }
    }
  }
}

static void putPrevInteractingInTable(int atom,struct parameters *par,struct systemPos *pos,int *table){
  int listPosj;
  
  listPosj=prevInteractionHead[atom]; 
  while (prevInteractionList[listPosj]!=-1)
    table[prevInteractionList[listPosj++]]=1;
}


static void putInteractingInPrevList(int atom,struct parameters *par,struct systemPos *pos,struct bondList *blist,  int *posInPrevList){
  int j,bondListPosj,ngbrListPosj;
  double drx_ij,dry_ij,drz_ij;
  double rij2;
  double r02;
  
  /* r02 is the distance at which ngbr atoms are included, shoud this include skin ? probably not! */
  /*  r02=POW2(sqrt(wwwPar.repR02)+par->skin); */
  r02=wwwPar.repR02;
  
  
  /* then the atoms it is bonded to */
  bondListPosj=blist->head[atom]; /*where the bonds start*/
  while(blist->list[bondListPosj]!=-1) /*get bonds of i, stops when the bonds od atom i+1 comes*/
    prevInteractionList[(*posInPrevList)++]=blist->list[bondListPosj++];
  
  
  /* and last the atoms it is ngbr with */
  if(wwwPar.includeRepPot){
    ngbrListPosj=nd.head[atom]; /*where the ngbrs start*/
    while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=nd.list[ngbrListPosj++]; /*j is the ngbr*/
	
      drx_ij=pos->x[j]-pos->x[atom]; 
      dry_ij=pos->y[j]-pos->y[atom]; 
      drz_ij=pos->z[j]-pos->z[atom];	
      PERIODIC(drx_ij,dry_ij,drz_ij);
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      if( rij2 <r02 )
		  prevInteractionList[(*posInPrevList)++]=j;
	
    }
  }
  prevInteractionList[(*posInPrevList)++]=-1;
}

/*
  Description:

  updateParts takes in two lists, part and partSkin and updates these so that the size of them are increased. This is done by creating a new part from the old part by including all the atoms that each atom in the old part is interacting with through the potential if the atom in the old part has moved by more than udateTreshold compared to the position it had prior to the starting of the current minimization (these are stored in pos.prevxyz). Partskin is then updated to be a list containing all atoms in the new part and all the atoms they are interacting with

  Input: 
  par             The parameter struct
  pos             The position struct
  part            An array containing atomsInPart atoms. 
  partSkin        An array containing atomsInPartSkin atoms
  atomsInPart     See part
  atomsInPartSkin See partSkin
  includePrevInteractionsInSkin    this include sin partskin also all the atoms that part atoms previously had som interactions with (at beginning of current minimization) 
  updatePart      if 0 then it does not update part (it always upates partSkin)

  Output:
  none

  Modified parameters:
  part
  partSkin,
  tomsInPart     

    atomsInPartSkin 
    */

static void updateParts(struct parameters *par,struct systemPos pos,int *part,int *partSkin, int *atomsInPart,int *atomsInPartSkin,int includePrevInteractionsInSkin,int updatePart){
  static int firsttime=1;
  int i,ip,j;
  struct bondList *blist;
  static int *table;
  double drx,dry,drz;
  double r2;
  
  
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  
  if(firsttime){
    table=malloc(pos.nAtoms*sizeof(int)); /* in this table we put a 1 if it should be in partSkin (the new one) */
    memset(table,0,pos.nAtoms*sizeof(int)); /*after this we have to make sure that table is initialize to 0 once we have used it since it is not done except the firsttime*/
    firsttime=0;
  }
  getBondList(&blist);
  

  /* first new part has to be created */
  
  if(updatePart){	    
    for(ip=0;ip<(*atomsInPart);ip++)
		table[part[ip]]=1;
	 for(ip=0;ip<(*atomsInPart);ip++)
		putInteractingInTable(part[ip],par,&pos,blist,table);/* lets  include its neighbours*/
  

	 
    /* now write this to the part list and put the table to zero again */
    j=0;
    for(i=0;i<pos.nAtoms;i++)
      if(table[i]){
		  part[j++]=i;
		  table[i]=0; /*make sure table only contains 0:s.. */
      }
    (*atomsInPart)=j;
  }

  /* then new partSkin has to be created */

  /* put all  interacting particles */
  for(ip=0;ip<(*atomsInPart);ip++)
		table[part[ip]]=1;

  for(ip=0;ip<(*atomsInPart);ip++)
    putInteractingInTable(part[ip],par,&pos,blist,table);
  
  /* this include sin partskin also all the atoms that part atoms previously had som interactions with (at beginning of current minimization) */
  if(includePrevInteractionsInSkin)	  
    for(ip=0;ip<(*atomsInPart);ip++)
      putPrevInteractingInTable(part[ip],par,&pos,table);
  
  /* now write this to the partSkin table */
  j=0;
  for(i=0;i<pos.nAtoms;i++)
    if(table[i]){
      partSkin[j++]=i;
      table[i]=0; /*make sure table only contains 0:s.. */
    }
  (*atomsInPartSkin)=j;
}

double ringPenaltyEnergy(struct parameters *par,struct systemPos pos){
  int path[10];
  int i,jlpos,klpos;
  struct bondList *blist;
  double pot=0;
  int pathlength;
  int numOfNgbrs;
  int listStart,listEnd;
  int n;
  getBondList(&blist);

  for(i=0;i<pos.nAtoms;i++) 
    if(pos.aType[i]==SI_ATYPE){
      /* firstl loop through all angles at SiAtomss */
      numOfNgbrs=numBonds(i);
      listStart=blist->head[i];
      listEnd=listStart+numOfNgbrs;
      
      for(jlpos=listStart;jlpos<listEnd;jlpos++)
		  for(klpos=jlpos+1;klpos<listEnd;klpos++){
			 n=getShortSilicaRing(blist,pos.aType,blist->list[jlpos],i,blist->list[klpos]);
			 if(n==3){
				pot+=0.3333333333*wwwPar.ringPenalty3;
			 }
			 else if(n==2)
				pot+=0.5*wwwPar.ringPenalty2;
	    
		  }
    }
  return pot;
}






void reinitWwwAfterBondChange(struct parameters *par,struct systemPos pos){
  if(wwwPar.includeRingPot)
    wwwPar.ringPenaltyEnergy=ringPenaltyEnergy(par,pos);
  else
    wwwPar.ringPenaltyEnergy=0.0;
}


/*
  reinitWwwPotAfterStep can be run if a step which is not accepted has been done, it simply copies a old version of the ngbrlist on top of the current one so that the next minimization can start from a correct situation 
  Also updates	prevInteractionHead and prevInteractionList
*/

void reinitWwwPotAfterStep(struct parameters *par,struct systemPos pos,int stepAccepted){
  int i,listPos;
  struct bondList *blist;
  getBondList(&blist);
  
  if(stepAccepted){
#ifdef PARALLEL
    updateNgbrs(&nd,pos,par);
    nonNnNgbrs(&nd,pos.nAtoms);     
#endif
    copyNgbrs(&ndLastStep,&nd,pos.nAtoms);
    listPos=0;
    for(i=0;i<pos.nAtoms;i++){
      prevInteractionHead[i]=listPos;
      putInteractingInPrevList(i,par,&pos,blist,&listPos);
    }

#ifdef PARALLEL  
    wwwPar.ringPenaltyEnergy=ringPenaltyEnergy(par,pos); /*if paralell, make sure this value is correct (cant use  wwwPar.ringPenaltyEnergy) */
#endif   
    prevRingPenaltyEnergy=wwwPar.ringPenaltyEnergy; 
  }
  else{
    copyNgbrs(&nd,&ndLastStep,pos.nAtoms); /* copy back to ngbr previous ngbr list */
    wwwPar.ringPenaltyEnergy=prevRingPenaltyEnergy;
  }
  
}












