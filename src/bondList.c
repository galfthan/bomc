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
#include "shared.h"
#include "mcSubs.h"
#include "bondList.h"
#include "randFunc.h"

/* this module keeps the bondlist. */

static int checkBondValidity(struct systemPos *pos, int i,int j);

static struct bondList blist;	  /* this is the bondlist that the rest of the prgoram uses */
static int *tempList;

  
void	getBondInfo(int **head,int **list){
  *head=blist.head;
  *list=blist.list;
}

void getBondList(struct bondList  **lblist){ /* used to get a pointer to the bondlist */
  *lblist=&blist;
}


/*
  this creates an empty bondlist
*/

void initBondList(struct systemPos pos){
  int i,j;
  
  
  tempList=malloc(sizeof(int)*pos.nAtoms); /* this is a temporary list used in this module by some functions */
  blist.nAtoms=pos.nAtoms;
  blist.nMax=5;
  blist.arraySize=2*pos.nAtoms+pos.nAtoms*blist.nMax;
  
  /*these are allocated as one big list after each other, this is so that we can easily transfer it from one
    processor to another*/
  
  blist.head=malloc(sizeof(int)* blist.arraySize);
  blist.danglBonds=blist.head+pos.nAtoms;
  blist.list=blist.head+2*pos.nAtoms;
  
  for(i=0;i<blist.nAtoms;i++){
    if(pos.aType[i]==SI_ATYPE)
      blist.danglBonds[i]=4;
    else  if(pos.aType[i]==O_ATYPE)
      blist.danglBonds[i]=2;
	 
    blist.head[i]=blist.nMax*i;
    for(j=0;j<blist.nMax;j++)
      blist.list[blist.head[i]+j]=-1;

  }
  

}
void deleteLongBonds(struct systemPos pos,struct parameters *par,double rmax){
  int i,ai,aj;
  double rij2,drx_ij,dry_ij,drz_ij;
  double rmax2;
  rmax2=rmax*rmax;
  for(ai=0;ai< pos.nAtoms;ai++){ /*loop over atoms i*/
    i=blist.head[ai];
    while((aj=blist.list[i++])!=-1){
      drx_ij=pos.x[aj]-pos.x[ai]; 
      dry_ij=pos.y[aj]-pos.y[ai]; 
      drz_ij=pos.z[aj]-pos.z[ai];	
      PERIODIC(drx_ij,dry_ij,drz_ij);
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      if(rij2>rmax2)
	removeBond(ai,aj);
    }
  }
}
/* onnectAtoms connects atoms with bonds. It does it by cconnectin all atoms with a distance less than rcut between them.  */
  void connectAtoms(struct systemPos pos,struct parameters *par,double rcut){
  int i,j,k,bondListPosj,bondListPosk;
  double rij2,drx_ij,dry_ij,drz_ij;
  int appendBond;
  double rik2,drx_ik,dry_ik,drz_ik;
  int totNumOfBonds=0;
  int numOfBonds=0;

  bondListPosj=0;
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    numOfBonds=0;
    
    blist.head[i]=bondListPosj;
    for(j=0;j< pos.nAtoms;j++) if(i!=j){ /*loop over atoms j*/
      drx_ij=pos.x[j]-pos.x[i]; 
      dry_ij=pos.y[j]-pos.y[i]; 
      drz_ij=pos.z[j]-pos.z[i];	
      PERIODIC(drx_ij,dry_ij,drz_ij);
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      
      if( !(pos.aType[i]==O_ATYPE && pos.aType[j]==O_ATYPE) && rij2<POW2(rcut)){ /*if not a O-O bond, accept it and distance smaller than rcut (in ngbrlist we have rcut+skin)*/
	
	/* first we ccheck that there are no bonds at atom i which have a low angle to ij, if there are we choose the shorter one of these */
	appendBond=1;
	bondListPosk=blist.head[i];
	while(bondListPosk<bondListPosj){
	  k=blist.list[bondListPosk];
	  drx_ik=pos.x[k]-pos.x[i]; 
	  dry_ik=pos.y[k]-pos.y[i]; 
	  drz_ik=pos.z[k]-pos.z[i];	
	  PERIODIC(drx_ik,dry_ik,drz_ik);
	  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
			
	  if((drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/sqrt(rik2*rij2)>0.95){ /* small angle between bonds ik and ij */
	    appendBond=0;	  /* here bond ij will not be appended last but will wither replace ik or not be inserted at all */
	    if(rik2>rij2)	  /* replace ik by ij */
	      blist.list[bondListPosk]=j;
	    break;
	  }
	  bondListPosk++;
	}
	if(appendBond){		  /* if there was no bond at low angle accept this new bond and put it to end of */
	  blist.list[bondListPosj++]=j;
	  numOfBonds++;
	  totNumOfBonds++;
	}
		  
			 
      } 
    }
	 
	 
    for(j=numOfBonds;j<blist.nMax;j++)	  /* we will put endmarkers so that all atoms have place for 4 bonds */
      blist.list[bondListPosj++]=-1;/*put endmarker*/
	 
  }
  
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    if(pos.aType[i]==SI_ATYPE)
      blist.danglBonds[i]=4-numBonds(i);
    else
      blist.danglBonds[i]=2-numBonds(i);
  }
  
  printf("connectatoms ready, found %d bonds %g per atom\n",totNumOfBonds,(double)totNumOfBonds/pos.nAtoms);
  
}
/*this functoion makes sure that no first or second ngbrs in bondlist are included in nngbrlist*/
void nonNnNgbrs(struct ngbrData *nd,int nAtoms){
  int i,j,jb,kb,ngbrListPosj,accNgbrListPosj,bondListPosj,bondListPosk;
  
  for(i=0;i< nAtoms;i++){ /*loop over atoms i*/
    INITLISTPOS(nd->head,i,ngbrListPosj);
    accNgbrListPosj= ngbrListPosj; /*accNgbrListPosj shows where we are in the ngbr list when we write the accepted atoms*/
    LOOPNGBR(nd->list,j,ngbrListPosj){
      /*lets check if the ngbr j is in the bondlist as a first or second ngbr, if not then we will calculate the energy of this bond*/
      
      INITLISTPOS(blist.head,i,bondListPosj); /*where the ngbrs start*/
      LOOPNGBR(blist.list,jb,bondListPosj){
	if(j==jb) /*j is the ngbr not accepted*/
	  goto notAccepted;
	
	INITLISTPOS(blist.head,jb,bondListPosk); /*here we start to check second ngbr*/
	LOOPNGBR(blist.list,kb,bondListPosk)
	  if(j==kb) /*if j is second ngbr not accepted*/
	    goto notAccepted;
      }
      nd->list[accNgbrListPosj++]=j;
      
    notAccepted:;
    }
    for(j=accNgbrListPosj;j<ngbrListPosj;j++)
      nd->list[j]=-1; /*put the endmarker for this list*/
    
  }
  
}
/*this functoion makes sure that no first or second ngbrs in bondlist are included in nngbrlist*/
void nonNnNgbrsPart(struct ngbrData *nd,int nAtoms,int *partSkin,int atomsInPartSkin){
  int i,ips,j,jb,kb,ngbrListPosj,accNgbrListPosj,bondListPosj,bondListPosk;
 
  //  for(i=0;i< nAtoms;i++){ /*loop over atoms i*/
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    INITLISTPOS(nd->head,i,ngbrListPosj);
    accNgbrListPosj= ngbrListPosj; /*accNgbrListPosj shows where we are in the ngbr list when we write the accepted atoms*/
    LOOPNGBR(nd->list,j,ngbrListPosj){
      /*lets check if the ngbr j is in the bondlist as a first or second ngbr, if not then we will calculate the energy of this bond*/
      
      INITLISTPOS(blist.head,i,bondListPosj); /*where the ngbrs start*/
      LOOPNGBR(blist.list,jb,bondListPosj){
	if(j==jb) /*j is the ngbr not accepted*/
	  goto notAccepted;
		  
	INITLISTPOS(blist.head,jb,bondListPosk); /*here we start to check second ngbr*/
	LOOPNGBR(blist.list,kb,bondListPosk)
	  if(j==kb) /*if j is second ngbr not accepted*/
	    goto notAccepted;
      }
      nd->list[accNgbrListPosj++]=j;
      
    notAccepted:;
    }
    for(j=accNgbrListPosj;j<ngbrListPosj;j++)
      nd->list[j]=-1; /*put the endmarker for this list*/
	 
  }
  
}

int numBonds(int ai){
  int i,ain=0;
  i=blist.head[ai];
  while(blist.list[i++]!=-1)
    ain++;
  return ain;
}

int getRandBondAtAtom(int ai){
  double r=randNum(2);
  return blist.list[blist.head[ai]+(int)(randNum(2)*numBonds(ai))];
}

int areBonded(int ai,int aj){
  int i,bondAi=0,bondAj=0;
  
  /*first check that they are bonded... */
  i=blist.head[ai];
  while(blist.list[i]!=-1)
    if(blist.list[i++]==aj) bondAi++;

  i=blist.head[aj];
  while(blist.list[i]!=-1)
    if(blist.list[i++]==ai) bondAj++;
  
  if(bondAi==1 &&bondAj==1 ) 
    return 1;						  /* Everything ok */
  else
    return 0;

}
/*!
  \brief Checks if ai and aj have a common ngbr to which theya re bonded.
	
*/
int areBonded2(int ai,int aj){
  int i;
  i=blist.head[ai];
  while(blist.list[i]!=-1){
    if(areBonded(aj,(blist.list[i])))
      return 1;
    i++;
  }
  return 0;							  /* not 2bonded */
}

/* a1 and a2 are bonded. This function will traverse starting from a1->a2 until it finds a atom with more than two bonds. This is then returned. Thus  a oxygen between two silicons can be jumped over automaitcally  */
int getNextNode(int a1, int a2){
  int listPos;
  int a3;
  
  while(numBonds(a2)==2){		  /* loop until we get a atom which has more than two bonds */
    listPos=blist.head[a2];  
    while((a3=blist.list[listPos++])!=-1){	  /* loop over bonds of a2 */
      if(a3!=a1){					  /* if bond not backwards then this is the next bond we move to */
	a1=a2;
	a2=a3;
	break;						  /* start again for this new a2 and see if this is ok, has more than two bonds */
      }
    }
    
  }
  return a2;
}




/*get a andom atom with dangling bonds*/
int getRandAtomWdb(int minDanglBonds){ 
  
  int dbAtoms,i;
  
  dbAtoms=0;
  for(i=0;i<blist.nAtoms;i++)
    if(blist.danglBonds[i]>=minDanglBonds)
      tempList[dbAtoms++]=i;
  
  if(dbAtoms==0) /*no atoms with at least minDanglBonds dangling bonds*/
    return -1;
  
  return tempList[(int)(dbAtoms*randNum(2))];
}


/*get a andom atom with bonds*/
int getRandAtomWb(int minBonds){ 
  int ai,i;
  int bAtoms=0;
	
  /* firt calculate atom with algorithm that is good if there are lots of atoms with bonds */
  for(i=0;i<100;i++){
    ai=(int)(blist.nAtoms*randNum(2));
    if(numBonds(ai)>=minBonds)
      return ai;					  /* found accepted atom */
  }

  /* if we came here there are few or none atoms with bonds, lets use another algorithm */
  
  for(i=0;i<blist.nAtoms;i++)
    if(numBonds(ai)>=minBonds)
      tempList[bAtoms++]=i;
  
  if(bAtoms==0) /*no atoms with at least minBonds  bonds*/
    return -1;
  
  return tempList[(int)(bAtoms*randNum(2))];
}


int createBond(struct systemPos *pos,int ai,int  aj){
  if(checkBondValidity(pos,ai,aj))
    return 1; /* bond cant be created */ 

  blist.list[blist.head[ai]+numBonds(ai)]=aj; /*add bonds last in list */
  blist.list[blist.head[aj]+numBonds(aj)]=ai;
  blist.danglBonds[ai]--; /*both lost one dangling bonds */
  blist.danglBonds[aj]--;
  return 0;
}


int removeBond(int ai, int aj){
  int i,delta;
   
  if(!areBonded(ai,aj))
    return 1;
  
  /*first we will take away aj from ai:s list*/
  i=blist.head[ai];
  delta=0;
  while(blist.list[i]!=-1){
    if(blist.list[i]==aj) 
      delta++; /*if at aj then jump over this*/
    blist.list[i]=blist.list[delta+i];
    i++;
  }
  blist.list[i]=-1; /*change last one to a -1 stopper*/

  /*then we will take away ai from aj:s list*/
  i=blist.head[aj];
  delta=0;
  while(blist.list[i]!=-1){
    if(blist.list[i]==ai) 
      delta++; /*if at aj then jump over this*/
    blist.list[i]=blist.list[i+delta];
    i++;
  }
  blist.list[i]=-1; /*change last one to a -1 stopper*/
  
  blist.danglBonds[ai]++; /*both got one dangling bonds */
  blist.danglBonds[aj]++; 
  return 0;
}

/*!
  checks if a bond can be created between atoms i and j. It can be created if ia nd j are not both oxygen atoms annd are not already bonded and have no common ngbrs
*/
static int checkBondValidity(struct systemPos *pos, int i,int j){
  int ib,jb,listPosi,listPosj;
  /* cant be itself obvioulsy */
  if(i==j)
    return 1;

  if(blist.danglBonds[i]==0 || blist.danglBonds[j]==0 )
    return 2; /* bond cant be created unless they both have dangling bonds */ 
  
  /* check that no O-O bond is created */
  if(pos->aType[i]==O_ATYPE && pos->aType[j]==O_ATYPE)
    return 3;
  
  /* check that not already bonded */
  if(areBonded(i,j))
    return 4;  /*no double bonds are permitted*/

  /* check that they have no common ngbrs, we do not want to create systems with three atoms in a ring 
     and also check that no ngbrs are connected =>no 4 rings.. in silica, in silicon they are allowed!! (should make it easier for system to find minimum? 
  */
  listPosi=blist.head[i];  
  while((ib=blist.list[listPosi++])!=-1){ /* loop over ai:s ngbrs */
    listPosj=blist.head[j];  
    while((jb=blist.list[listPosj++])!=-1){ /* loop over ai:s ngbrs */
      if(jb==ib) 
	return 5;
      if(areBonded(ib,jb) && (pos->aType[i]==O_ATYPE || pos->aType[j]==O_ATYPE))
	return 6;
	
    }
  }
  

  return 0;
}
