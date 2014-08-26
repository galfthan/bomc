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
#include <string.h>

#include "shared.h"
#include "mcSubs.h"
#include "miscSubs.h"
#include "bondList.h"
#include "potentialWww.h"
#include "conjgrad.h"
#include "rings.h"

static int recursiveShortestRingCheck(struct bondList *blist, int *path,int atom,int pathDepth,int realDepth,int maxDepth); /* used by getRing, better not touch it otherwise */

static int destroyRing(struct systemPos pos,struct parameters *par, int ringlength,int *path, int pathlength);

/*!
  This destroys rings of length length by doing a bondswitch that destroys the ring while increasing the energy with as little as possible
  

  int destroyRings(int length){}
*/

/*! 
  This function creates a histogram of rings. 
  It is normalized by dividing with N, number of silicon atoms. of the rings to stdout
  stattype:1=For each angle count the shortest ring =>6N rings
  Statistics not divided with ringsize=>same ring counted multiple time (like in PRB 47,3053)
  2=Only shortest ring=> 6N rings.
  3=irreducible rings (all shortest rings are counted) 
  4=reducible rings   (all rings are counted...)       
			  

*/


int destroyRings(struct systemPos pos,struct parameters *par,int length){
  int *path;
  int destroyedRings=0;
  int failedDestruction=0;
  int pathlength;
  int i,j,k;
  int listposj,listposk;
  struct bondList *blist;
  
  getBondList(&blist);
  path=malloc(sizeof(int)*2*(length+1));
  
  

  for(i=0;i<pos.nAtoms;i++) if(numBonds(i)>2){ /* firstl loop through all angles at atoms with more than two bonds */
  startAgain:	
    listposj=blist->head[i];
    while((j=blist->list[listposj++])!=-1){
      listposk=listposj;
      while((k=blist->list[listposk++])!=-1){
       	pathlength=getRing(path,j,i,k,length);
	if(pathlength>0){ /* if  there is a ring  with this length this is true, now we should delete it!*/
	  int error=destroyRing(pos,par,length,path,pathlength);
	  if(!error){
	    destroyedRings++;
	    goto startAgain;	  /* the problem here is that destroyring changes the lists we are looping through, 
				     if we are unlucky (and sometimes we are) j might not be bonded to i anymore, better to start again from the top
				     sorry about the goto..
				  */
	  }
	  else{
	    failedDestruction++;
	  }
	}
      }
    }
  }

	 

  printf("DestroyedRings %d failed %d \n",destroyedRings,failedDestruction);
  free(path);
  return failedDestruction;
}


static void getPossibleSwitchesForDestruction(struct systemPos pos,struct parameters *par,int ringlength,int *path, int pathlength,int *switchList, int *possibleSwitches){
  int i,j;
  int tempPath[100];
  int ai,aj,aib,ajb;
  int posInBlist;
  int posInajBlist;
  int error;
  int prevInRing,nextInRing;
  struct bondList *blist;
  getBondList(&blist);
  (*possibleSwitches)=0;

  for(i=0;i<pathlength;i++){
    if(numBonds(path[i])>2){  /* this is a atom me might want to switch */
      ai=path[i];					  /* here is ai the first to participate in switch */
		
      if(i==pathlength-1) nextInRing=0; /* calculate taking into account the next and previous atom on ring */
      else                nextInRing=i+1;
		
      if(i==0) prevInRing=pathlength-1;
      else     prevInRing=i-1;
		

      for(j=0;j<2;j++){			  /* loop over aib (in a ugly way...) */
	if(j==0) aib=path[prevInRing];
	else     aib=path[nextInRing];
		  
	posInBlist=blist->head[ai];
	while((aj=blist->list[posInBlist++])!=-1) /* loop over ngbrs of ai  to choose aj*/
	  if(aj!=path[prevInRing] && aj!=path[nextInRing]) { /* make sure it is  not on ring */
	    aj=getNextNode(ai,aj); /* now we jump over a possible oxygen atom, makes sure that aj has more than two bonds */
				
	    posInajBlist=blist->head[aj];
	    while((ajb=blist->list[posInajBlist++])!=-1) /* loop over ngbrs of aj =ajb*/
	      if(!areBonded(ai,ajb)){ /* make sure ajb and ai not ngbrs */
		error=makeSpecificSwitch(&pos,ai,aj,aib,ajb); /* try to make theswitch so that we can check if it is  acceptable */
					 
		/* check that switch was successful and that these new bonds are not part of a too small ring (it happens in some special cases..) */
		if(!error && 
		   getBondMinRing(tempPath,ai,ajb,ringlength)==-1 &&
		   getBondMinRing(tempPath,aj,aib,ringlength)==-1){
						
		  switchList[*possibleSwitches*4+0]=ai; /* write down the switch*/
		  switchList[*possibleSwitches*4+1]=aib;
		  switchList[*possibleSwitches*4+2]=aj;
		  switchList[*possibleSwitches*4+3]=ajb;
		  (*possibleSwitches)++;
						
		}
		restoreState(pos); /* restore... , we do not want to make the switch just yet */
					 
	      }
	  }
      }
    }
  }
}

static int destroyRing(struct systemPos pos,struct parameters *par, int ringlength,int *path, int pathlength){
  int error;
  int i,j;
  int *switchList;
  int switches;
  double potVal;
  double minPot=1e30;
  int minSwitch=-1;
  static int *part;
  static int firsttime=1;
  int atomsInPart=4;
  
  if( firsttime){
    part=malloc(pos.nAtoms*sizeof(int));;
    firsttime=0;
  }
 
  switchList=malloc(sizeof(int)*ringlength*2*3*4*2); /* each atom on ring has 4-2=2 pssible bonds.
							these have three additional bonds gives length*2*3 switches 
							which each require 4 int:s to store. 
							Last 2 a precaution if this is somehow wrong (hate signal 11 :) */

 
  
  getPossibleSwitchesForDestruction(pos,par,ringlength,path,pathlength,switchList,&switches);
  printf("Removing ring, -");
  for(i=0;i<pathlength;i++)
    printf("%d-",path[i]);
  printf("  (%d possible switches) ...",switches);

  fflush(stdout);
  if(switches==0)
    return 10;
  
  for(i=0;i<switches;i++){
    error=makeSpecificSwitch(&pos,switchList[i*4],switchList[i*4+2],switchList[i*4+1],switchList[i*4+3]); /* make the switch itself */
    if(!error){
      for(j=0;j<4;j++)  part[j]=switchList[i*4+j];
      potVal=initialPartCg( pos,par,part,&atomsInPart,minPot,partwwwPot,wwwPot); /* here we can use optimization by giving minPot as treshold value! */
      if(potVal<minPot){
	minPot=potVal;
	minSwitch=i;
      }
    }
    restoreState(pos);
	 
  }
  
  if(minSwitch==-1)				  /* no switch was performed... */
    return 11;
  
  makeSpecificSwitch(&pos,switchList[minSwitch*4],switchList[minSwitch*4+2],switchList[minSwitch*4+1],switchList[minSwitch*4+3]);
  potVal=fullCg(pos,par,0,0,1e10,wwwPot); /*optimization not allowed, accPotVal has a arbitrary value */
  backupState(pos);
  
  printf("Ready. pot: %gev ",potVal/pos.nAtoms*TOEV);
  printf("Switch: %d-%d--%d-%d \n",switchList[minSwitch*4+1],switchList[minSwitch*4],switchList[minSwitch*4+2],switchList[minSwitch*4+3]); 
  return 0;
  
}



/*! calcualtes minimum ringsize of al-ar which is less than maxLength (that is the min ring this vertex is part of) */
int getBondMinRing(int *path,int a1,int a2,int maxLength){
  int i,a3,listPos;
  int pathlength;
  struct bondList *blist;
  getBondList(&blist); 
  
  if(!areBonded(a1,a2))
    return -3;
	  
  
  for(i=2;i<=maxLength;i++){
    /* first the rings centered at a1 */
    listPos=blist->head[a1];
    while((a3=blist->list[listPos++])!=-1) 
      if(a3!=a2 && (pathlength=getRing(path,a3,a1,a2,i))!=-1)
	return pathlength;
    /* then similarly the rings cetered at a2 */
    listPos=blist->head[a2];
    while((a3=blist->list[listPos++])!=-1) 
      if(a3!=a1 &&  (pathlength=getRing(path,a3,a2,a1,i))!=-1)
	return pathlength;
  }
  
  return -1;			  /* did not find */
}

/*! calcualtes minimum ringsize of al-ac-al which is less than maxLength   (that is the min ring this angle is part of)*/
int getMinRing(int *path,int ar,int ac,int al,int maxLength){
  int i;
  int pathlength;
  if(!areBonded(ar,ac)|| !areBonded(ac,al))
    return -2;
  
  for(i=2;i<=maxLength;i++)
    if((pathlength=getRing(path,ar,ac,al,i))!=-1)
      return pathlength;
  
  return -1;						  /* did not find */
  
}



/*! this function gets ring of length length starting at al-ac-ar and puts it in path, 
  returns -1 if no path was found, otherwise returns the amount of atoms in ring (might be other than length!) */


int getRing(int *path,int ar, int ac,int al,int length){
  struct bondList *blist;
  int realDepth=0;
  if(!areBonded(ar,ac)|| !areBonded(ac,al))
    return -2;
  
  
  getBondList(&blist);
  path[0]=ar;
  path[1]=ac;
  if(numBonds(ar)>2) realDepth++;
  if(numBonds(ac)>2) realDepth++;
  return recursiveShortestRingCheck(blist,path,al,2,realDepth,length);
}


int getShortSilicaRing(struct bondList *blist,int *atype,int ar, int ac,int al){
  int arNgbrs[4];
  int i,j,k;
  int numOfarNgbrs;
  int listStart,listEnd;
  int listpos,ngbrlistpos;
 
  if(atype[ar]!=O_ATYPE ||atype[al]!=O_ATYPE  ||atype[ac]!=SI_ATYPE)
    return 0;

  ar=getNextNode(ac,ar);
  al=getNextNode(ac,al);
  
  
  if(ar==al)
    return 2; /*this is a silica 2 ring*/

  numOfarNgbrs=numBonds(ar);
 
  listStart=blist->head[ar];
  listEnd=listStart+numOfarNgbrs;
  
  ngbrlistpos=0;
  for(listpos=listStart;listpos<listEnd;listpos++)
    if(atype[blist->list[listpos]]==O_ATYPE)
      arNgbrs[ngbrlistpos++]=getNextNode(ar,blist->list[listpos]);
  
  numOfarNgbrs=ngbrlistpos;
  
   
  for(j=0;j<numOfarNgbrs;j++)
    if(al==arNgbrs[j])
      return 3; /*is in a silica 3 ring*/
  
  return 0;
}

/*! 
  this function returns the number of atoms in the  ring of length maxDepth  that contains atom and its two ngbrs ngbr1 and ngbr2, return -1 if there are no riings with this length. 

  pathDepth is the amount of atoms in the path. RealDepth is the amount o atoms in path that have more than two ngbrs. If atom has two ngbrs then it does not participate in the 

  returns the pathDepth of found path
*/
static int recursiveShortestRingCheck(struct bondList *blist, int *path,int atom,int pathDepth,int realDepth,int maxDepth){

  int ngbrlpos;
  int ngbr;
  int n,i;
  int accepted;
  int numOfNgbrs;
  int listStart,listEnd;

  numOfNgbrs=numBonds(atom);
  if(numOfNgbrs>2) /* if more than two bonds then increase realDepth  */
    realDepth++;
  
  if(realDepth>maxDepth)				  /* too deep ,abort */
    return -1;
  
  path[pathDepth]=atom;			  /* add atom to path and increase pathlength */
  pathDepth++;						  
  
  
  listStart=blist->head[atom];
  listEnd=listStart+numOfNgbrs;
  
  
  for(ngbrlpos=listStart;ngbrlpos<listEnd;ngbrlpos++){  /* loop through ngbrs of atom and check the possible paths */
    ngbr=blist->list[ngbrlpos];
    if(realDepth==maxDepth && ngbr==path[0])	{	  /* we found it, the other end of the path!! , and this is the correct depth*/
      return pathDepth;
    } 
    
    /* here we first make sure that our path does not intersect itself */
    accepted=1;			
    for(i=pathDepth-1;i>0;i--) 
      if(ngbr==path[i]){
	  accepted=0;
	  break;
      }
    
    if(!accepted) continue;
    
    n=recursiveShortestRingCheck(blist,path,ngbr,pathDepth,realDepth,maxDepth); /* check next depth */
    if(n>0)
      return n; 	/* somebody after us found it!! */
	 
  }
  
  return -1;						  /* nah, did not ind anything, sorry */
}
