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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>
#include "shared.h"
#include "ngbrs.h"
#include "mcSubs.h"
#include "miscSubs.h"
#include "parallel.h"
#include "randFunc.h"
#include "bondList.h"
#include "initialize.h"
#include "potentialWww.h"
#include "fileio.h"

#define MAXBONDCHANGE 100
/*in this struct trial changes are stored*/
struct bondChanges{
  int numOf; 
  int *addr[MAXBONDCHANGE]; 
  int value[MAXBONDCHANGE];
};


static struct bondChanges bondc;
static int siType=1,oType=0;



static void backupBonds(int ai){
  int listPos;
  struct bondList *blist;
  getBondList(&blist);

  listPos=blist->head[ai];
  do{ /*put all bond in backup list and also the last -1 stop value*/
    bondc.addr[bondc.numOf]=&(blist->list[listPos]);
    bondc.value[bondc.numOf]=blist->list[listPos];
    bondc.numOf++;
  } while(blist->list[listPos++]!=-1);
  
  bondc.addr[bondc.numOf]=&(blist->danglBonds[ai]);
  bondc.value[bondc.numOf]=blist->danglBonds[ai];
  bondc.numOf++;
}

static void restoreBonds(){
  int i;
  struct bondList *blist;
  getBondList(&blist);
  for(i=0;i<bondc.numOf;i++)
    *(bondc.addr[i])=bondc.value[i];
  bondc.numOf=0;
  
}

void backupState(struct systemPos pos){
  memcpy(pos.prevx,pos.x,sizeof(double)*pos.nAtoms*3); /*save old pos*/  
  memcpy(pos.prevEnPot,pos.enPot,sizeof(double)*pos.nAtoms); /*save old pot*/  
  bondc.numOf=0; /*we backup this stateas the one we can fallback to =>no bonds to restore*/
}


void restoreState(struct systemPos pos){
  memcpy(pos.x,pos.prevx,sizeof(double)*pos.nAtoms*3); /*restore old pos*/
  memcpy(pos.enPot,pos.prevEnPot,sizeof(double)*pos.nAtoms); /*restore old pot*/  
  restoreBonds(); 
}

void mcVolume(struct systemPos pos,struct parameters *par,int restore){
  double scale[3];
  static double  prevbox[3];
  int i,j;
  
  if(restore==0){
    for(i=0;i<3;i++)    
      scale[i]=1.0+(randNum(2)-0.5)*0.01; /*lets allow up to  of change for now*/
    

    for(i=0;i<3;i++)
      if(par->varVol[i]){
		  prevbox[i]=par->box[i];
		  par->box[i]*=scale[i];
      }
    
    setSystemSize(par);
    bondc.numOf=0;
    
    
    for(j=0;j<3;j++)				  /* scale atom positions */
      if(par->varVol[j]) 
		  for(i=j*pos.nAtoms;i<(j+1)*pos.nAtoms;i++)
			 pos.x[i]*=scale[j];
    
    
    
  }
  else{
    /*pos should have been taken care of already in main*/
    for(i=0;i<3;i++)
      if(par->varVol[i])
	par->box[i]=prevbox[i];
    
    setSystemSize(par);
  }
}



/* joins two dangling bonds */
int mcBondCreate(struct systemPos *pos,struct parameters *par, int *part,int *atomsInPart){
  int ai,aj;
  double drx,dry,drz,r2,r02=25.0;
  int error;
  int count;

  ai=getRandAtomWdb(1);
  if(ai==-1) /*if -1 then there are no aotms with dangling bonds*/
	 return 1;
  count=0;
  do{
    if(count++>1000) return 20;
    aj=getRandAtomWdb(1);
  } while( ai==aj );
  
  if(pos->aType[ai]==oType  &&  pos->aType[aj]==oType)
	 return 2;
  
  
  drx=pos->x[ai]-pos->x[aj]; 
  dry=pos->y[ai]-pos->y[aj]; 
  drz=pos->z[ai]-pos->z[aj];	
  PERIODIC(drx,dry,drz);
  r2=LENGTH2(drx,dry,drz);
  if(r2>r02)
	 return 3; /*too far away, probably wont be accepted anyway, it is faster to reject this step here instead of doin an energy minimization*/

  /* ai and aj have been accepted if we are at this point */

  /* first backup some old values */
  bondc.numOf=0;
  backupBonds(ai);
  backupBonds(aj);
  
  /*then create bond*/
  
  error=createBond(pos,ai,aj);
  if(error){
	 restoreBonds();
	 return 4;  /*some error happened, perhaps double bond ?*/
  }
  
  part[0]=ai;
  part[1]=aj;
  (*atomsInPart)=2;

  return 0;  
}


int mcBondBreak(struct systemPos *pos, int *part,int *atomsInPart){
  int error;
  int ai,aj;
  struct bondList *blist;
  getBondList(&blist);

  ai=getRandAtomWb(1); /*first atom chosen at random among the ones with bonds*/
  aj=getRandBondAtAtom(ai);
  
  /*take backup*/
  bondc.numOf=0;
  backupBonds(ai);
  backupBonds(aj);
  
  
  error=removeBond(ai,aj);
  if(error){
	 restoreBonds();
	 return 1;
  }
  
  part[0]=ai;
  part[1]=aj;
  (*atomsInPart)=2;
  
  return 0;							  /* successful */
}



/*! 
	in  mcBondDiffuse  a dangling bond is moved from one atom (ai) to another one (aj)

	aj---ajb         aj~   ajb
	 \        ==>     \   /
	  ai~               ai
 
where ~ is a dangling bond.
	
 */


int mcBondDiffuse(struct systemPos *pos, int *part,int *atomsInPart){
  int ai,aj,ajb;
  int error;
  int count;
  
  
  ai=getRandAtomWdb(1); /*first atom chosen at random among the ones with dangling bonds*/
  
  if(ai==-1 || numBonds(ai)==0) /* there are no atoms with dangling bonds or the number of bonds on ai is zero=>abort */
	 return 1;
		
  aj=getRandBondAtAtom(ai); /*aj is ngbr of ai, has to have 2 or more bonds*/
  if(numBonds(aj)<2)
	 return 2;

  count=0;
  do{
    if(count++>1000) return 20;
	ajb=getRandBondAtAtom(aj);
  } while(ajb==ai); /*ajb and ai have to be different*/
  
  /*take backup*/
  bondc.numOf=0;
  backupBonds(ai);
  backupBonds(aj);
  backupBonds(ajb);
 
  /*do moves*/
  
  error=removeBond(aj,ajb);
  error+=createBond(pos,ai,ajb);
  if(error){
	 restoreBonds(); 
	 return 3;
  }

  
  part[0]=ai;
  part[1]=aj;
  part[2]=ajb;
  (*atomsInPart)=3;
  
  return 0;							  /* successful */
}

/* A oxygen can with this diffuse between intersilicon sites*/
int mcOxygenDiffuse(struct systemPos *pos, int *part,int *atomsInPart){
  int oxAtom,siAtom0,siAtom1,siAtom2;
  int error;
  int count;

  /*get atom with at least two bonds*/
  oxAtom=getRandAtomWb(2);
  if(oxAtom==-1)
    return 1;						  /* not accepted, could not find such a atom */
  
  
  if(pos->aType[oxAtom]!=oType)	  /* if not oxygen atom */
    return 1; /* not successful, we by chance did not pick a oxygen atom.
		 might ths be made better ?*/
  
  siAtom1=getRandBondAtAtom(oxAtom);	  /* picks a Silicon atom.. */
  /*now lets pick the other atom the oxygen atom is connected to, ugly while hack*/
  count=0;
  do{
    if(count++>1000) return 20;
    siAtom0=getRandBondAtAtom(oxAtom);	  /* picks a Silicon atom.. */
  } while(siAtom0==siAtom1);
  
  siAtom2=getRandBondAtAtom(siAtom1);	  /* picks a Silicon atom.. */
  
  if(pos->aType[siAtom2]==oType)
    return 2; /*the second silicon atom is a oxygen atom not good...*/
  
  /*if we got this far everything should be clear */
  
  /*take backup*/
  bondc.numOf=0;
  backupBonds(oxAtom);
  backupBonds(siAtom0);
  backupBonds(siAtom1);
  backupBonds(siAtom2);
  
  /*do moves*/
  
  error=removeBond(oxAtom,siAtom0);
  error+=removeBond(oxAtom,siAtom1);
  error+=removeBond(siAtom1,siAtom2);

  error+=createBond(pos,oxAtom,siAtom2);
  error+=createBond(pos,oxAtom,siAtom1);
  error+=createBond(pos,siAtom0,siAtom1);
  if(error){
	 restoreBonds(); 
	 return 3;
  }

  
  part[0]=siAtom0;
  part[1]=siAtom1;
  part[2]=siAtom2;
  part[3]=oxAtom;
  (*atomsInPart)=4;
  /*one could perhpas at this point change the position of the oxygen atom so that it would be halfway between the silicon atoms, one shopuld check though that no backups are taken of the postitons after this.... dont remember exactly anymore... */

  return 0;							  /* successful */
}


/*! \brief Makes a bond switch trial move


 */


/*this is a modified switch where one can choose atoms which are not closest ngbrs*/
int mcBondSwitch(struct systemPos *pos, int *part,int *atomsInPart){
   int ai,aj,aib,ajb;
   int error;
 
  struct bondList *blist;
  getBondList(&blist);
  
 
  /*lets choose atoms ai and aj which are connected.
    These atoms will exchange  a ngbr bond */
  
  ai=getRandAtomWb(2);
  if(ai==-1)
    return 1;		          /* not accepted, could not find such a atom */
  if(pos->aType[ai]==oType)	  /* if oxygen atom */
    ai=getRandBondAtAtom(ai);	  /* picks a Silicon atom.. */
  
  aj=getRandBondAtAtom(ai);  /*aj is ngbr of ai*/
  aj=getNextNode(ai,aj);     /*jump over possible oxygen atoms*/
  
  /* choose atoms to switch */
  aib=getRandBondAtAtom(ai);
  ajb=getRandBondAtAtom(aj);
  
  /*check if these are valid */
  if(aib==ajb || areBonded(aib,aj) ||areBonded(ajb,ai))
    return 2;
  
  error=makeSpecificSwitch(pos,ai,aj,aib,ajb);
  if(error) 	 
    return error;
  
  /* we succeeded ! lets make the part that is used in minimization */
  part[0]=ai;
  part[1]=aj;
  part[2]=aib;
  part[3]=ajb;
  (*atomsInPart)=4;
  return 0;
}



int makeSpecificSwitch(struct systemPos *pos,int ai,int aj, int aib, int ajb){
  int error;
  int aiBond;
  int ajBond;
  /*take backup*/
  bondc.numOf=0;
  backupBonds(ai);
  backupBonds(aj);
  backupBonds(aib);
  backupBonds(ajb);
  /*now the bonds have been chosen, lets actually change them*/
  /*lets switch the bonds*/
  
  /* first we check actually to which atom ai is bonded effeectively, that is if ajb only has two ngbrs(is oxygen) it goes to the next one along bond */
  aiBond=getNextNode(ai,aib);
  ajBond=getNextNode(aj,ajb);
  
  /* remove bonds */
  error=removeBond(aj,ajb);
  error+=removeBond(ai,aib);
  
  if(error){
    restoreBonds(); 
    return 5;
  }
  /* create new bonds */
  error=createBond(pos,ai,ajb);
  error+=createBond(pos,aj,aib);
  if(error){
    restoreBonds(); 
    return 6;
  }
 
  /* 
     if aib,ajb has two ngbr and ai is now effectivly connected to the same atom as over ajb as previously over aib
      then the topology has not changed, no sense in making this move!! 
  */   
  if(numBonds(ajb)==2 &&  numBonds(aib)==2 && aiBond==getNextNode(ai,ajb) && ajBond==getNextNode(aj,aib)){
    restoreBonds(); 
    return 7;
  }
  return 0;
}


double getKT(int iteration,int nAtoms){
  static int firsttime=1;
  double kt;
  static int annealLength,startAnneal;
  static double  annealScale,endKt,startKt;

  if(firsttime){
	 int error;
	 double temp;
	 firsttime=0;
	 openParameterFile("main.par");
	 error=getParValue("startKt",&startKt,"%lf");
	 startKt*=EV;
	 error+=getParValue("endKt",&endKt,"%lf");
	 endKt*=EV;
	 error+=getParValue("annealScale",&annealScale,"%lf");

	 error+=getParValue("annealLength",&temp,"%lf");
	 annealLength=(int)(temp*nAtoms)+1;
	 error+=getParValue("startAnneal",&temp,"%lf");
	 startAnneal=(int)(temp*nAtoms);
	 closeParameterFile();
	 if(error>0){
		printf("could not read cooling schedule from main.par\n");
		exit(0);
	 }
  
  }
  if(iteration<startAnneal)
	 kt=startKt;
  else
	 kt=pow(annealScale,(double)((iteration-startAnneal)/annealLength))*startKt; /* pow is slow but this is only one calculation so... */
  
  
  if(kt<endKt)
	 return -1.0;					  /* give negative value if the annealing is finished */
  else 
	 return kt;
  

}

