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
#include <memory.h>
#include "shared.h"
#include "ngbrs.h"
#include "miscSubs.h"


/*this is a module that calculates the ngbrs of atoms. the ngbr information is stored in a struct ngbrData
  which is defined in the .h file. This structure must first be initialized with initiNgbrs(), after that 
  one can update the ngbr list with updateNgbrs(). updateNgbrs() checks if we need to update the 
  ngbrlist so its a good idea to call updateNgbrs every simulationstep.
  
  struct ngbrData
    int *list;           contains the ngbrlists for the atoms, each separate list is ended with a -1 to mark the end
    int *head;           contains the starting position of the ngbrlist on list
    double *prevx,*prevy contains the positions of the atoms at the previous step, used to check
			 ,*prevz		      if we need to update the ngbrlist
    int nMax;            maximum number of ngbrs, it is automaically increased if neede
    double rcs;          rcut +skin,defines the maximum distance of the ngbrs from the atom
    struct cellData *cd; the data for the cellmethod which is used to speed up the ngbr calculation
    int forceUpdate;     if 1, the ngbrs are updated without checking if it is neccessary
   
	 List contains the list of the ngbrs of a certain atom. The last element in a list is -1.
	 head contains the information where the list begins for a certain atom. Right now the lists are in order, previously they were noot..
	 
	 
*/


int checkNgbrUpdateAll(struct ngbrData *nd, struct systemPos pos,struct parameters *par){
  int i;
  double rx,ry,rz,r2;
  double maxHr2;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
 
  /*check if we need to update*/
  if( nd->forceUpdate ){
	 nd->forceUpdate=0;
	 return 1;
  }

  maxHr2=0.0; /*higher max*/

  for(i=0;i<pos.nAtoms;i++){
    rx=pos.x[i]-nd->prevx[i];
    ry=pos.y[i]-nd->prevy[i];
    rz=pos.z[i]-nd->prevz[i];
    
    r2=LENGTH2(rx,ry,rz);
    if(r2>minhBox2){
      PERIODIC(rx,ry,rz);
      r2=LENGTH2(rx,ry,rz);
    }
    
    
    if(r2>=maxHr2){
      maxHr2=r2;
      if(par->skin <= 2.0*sqrt(maxHr2)){
		  return 1;
      }
    }
  }
  
  return (par->skin <= 2*sqrt(maxHr2));
} 

int checkNgbrUpdatePart(struct ngbrData *nd, struct systemPos pos,struct parameters *par,int *partSkin,int atomsInPartSkin){
  int i,ips;
  double rx,ry,rz,r2;
  double maxHr2;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  

  /*check if we need to update*/
  if( nd->forceUpdate ){
	 nd->forceUpdate=0;
	 return 1;
  }

  maxHr2=0.0; /*higher max*/
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    rx=pos.x[i]-nd->prevx[i];
    ry=pos.y[i]-nd->prevy[i];
    rz=pos.z[i]-nd->prevz[i];
    
    r2=LENGTH2(rx,ry,rz);
    if(r2>minhBox2){
      PERIODIC(rx,ry,rz);
      r2=LENGTH2(rx,ry,rz);
    }
    
	 
    if(r2>=maxHr2){
      maxHr2=r2;
      if(par->skin <= 2*sqrt(maxHr2))
	return 1;
    }
  }
  
  return (par->skin <= 2*sqrt(maxHr2));
} 


void initNgbrs(struct ngbrData *nd, struct systemPos pos,struct parameters *par,double rcut){
  nd->nMax=5; /*will be increased if needed by updateNgbrs*/
  nd->rcs=rcut+par->skin;
  nd->forceUpdate=1; /*force a update the first time, it is turned of each time a forced update is performed*/
  nd->cd=malloc(sizeof(struct cellData));
  nd->cd->lCell=nd->rcs;
  cellInit(nd->cd,pos,par,1);
  
  nd->list=malloc(sizeof(int)*pos.nAtoms*nd->nMax);
  nd->head=malloc(sizeof(int)*pos.nAtoms); 
 
  nd->prevx=malloc(sizeof(double)*3*pos.nAtoms);
  nd->prevy=nd->prevx+pos.nAtoms;
  nd->prevz=nd->prevx+2*pos.nAtoms;
  nd->volume=par->volume;
   
  /*init some values for checking of next update*/
  
  if(nd->list==NULL ||  nd->head==NULL ){
    printf("ERROR in allocating memory in initNgbrs()\n");
    exit(0);
  }

  memcpy(nd->prevx,pos.x ,3*pos.nAtoms*sizeof(double));
  updateNgbrs(nd,pos,par);
}

void clearNgbrs(struct ngbrData *nd){
  clearCells(nd->cd);
  free(nd->list);
  free(nd->head);
  free(nd->prevx);
}
void copyNgbrs(struct ngbrData *dest,struct ngbrData *source,int nAtoms){
  if(dest->nMax<source->nMax) {
    dest->list=realloc(dest->list,(size_t)sizeof(int)*nAtoms*source->nMax);
    dest->nMax=source->nMax;
  }
  memcpy(dest->list,source->list,sizeof(int)*source->nMax*nAtoms);
  memcpy(dest->head,source->head,sizeof(int)*nAtoms);
  memcpy(dest->prevx,source->prevx,sizeof(double)*3*nAtoms);
  dest->rcs=source->rcs;
  dest->forceUpdate=source->forceUpdate;
  copyCells(dest->cd,source->cd,nAtoms);
}


void updateNgbrs(struct ngbrData *nd, struct systemPos pos,struct parameters *par){
  
  int n,i,j,k;
  int posInList;
  double rx,ry,rz,r2,rcs2;
 
  static int updates=0;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  
  /*We are now updating the lists */
  updates++;
  /*init some values for checking of next update*/
 
  nd->forceUpdate=0;
  
  rcs2=nd->rcs*nd->rcs;  /* == (rc+rskin)*(rc+rskin);*/
 
  enforcePeriodicity(pos,par);  /* make sure atoms are placed periodically */
  
  /* take copies of position so that we can calculate when we shoould update ngbrs the next time*/
  memcpy(nd->prevx,pos.x ,3*pos.nAtoms*sizeof(double));
  
  if(nd->volume!=par->volume){
	 nd->volume=par->volume;
	 cellInit(nd->cd,pos,par,0);
  }
  /*update cells first*/
  cellUpdate(nd->cd,pos,par);
  
  posInList=0;
  for(i=0;i<pos.nAtoms;i++){ 

	 n=nd->cd->cellList[i]; /*n is the cell of atom i*/
	 nd->head[i]=posInList; /*save position where the ngbrs of atom i begins*/
	 for(k=0;k<27;k++) /*loop over ngbrcells*/
		if(nd->cd->ngbrCells[n*27+k]>=0) /*if it is ok*/
		  for(j=nd->cd->head[nd->cd->ngbrCells[n*27+k]]; j >=0 ; j = nd->cd->list[j]){ /*loops over atoms in ngbrcell*/
			 
			 rx=pos.x[i]-pos.x[j];
			 ry=pos.y[i]-pos.y[j];
			 rz=pos.z[i]-pos.z[j];
					 
			 r2=LENGTH2(rx,ry,rz);
			 
			 if(r2>minhBox2){
				PERIODIC(rx,ry,rz);
				r2=LENGTH2(rx,ry,rz);
			 }
			 
			 if( r2 <= rcs2 && j!=i){
				nd->list[posInList  ]=j;
				posInList++;
			 }
				
			 /*if the ngbrs dont fit in list increase list size (-1 so that endmarker will always fit)*/
			 if(posInList>=pos.nAtoms*nd->nMax-10){
				nd->nMax=nd->nMax*1.5+1;
				nd->list=realloc(nd->list,(size_t)sizeof(int)*pos.nAtoms*nd->nMax);
			 }
		  }
	 /*leave marker for end of list*/
	 nd->list[posInList++]=-1;
	 
  } 


  
  if(DEBUG && updates%1000==0){
	 printf("ngbrlist updates %d\n",updates);
	
  }
}

/*
  This updates the ngbrlists of atoms on partSkin. One importatn note is that the aroms outside this list keeps their original 
  lists. This means that atoms outside lists can have wrong ngbrlist IF they are near enough to atoms which have moved. This can be avoided by making sure that partSkin includes ALL atoms that are close enough to the atoms that have moved..
  

*/
void updateNgbrsPart(struct ngbrData *nd, struct systemPos pos,struct parameters *par,int *partSkin,int atomsInPartSkin){
  
  int n,i,ips,j,k;
  int posInList,ndlistPos;
  
  double rx,ry,rz,r2,rcs2;
  static int firsttime=1;
  static int tempListnMax;
  static int *tempList;
  static int *tempHead;
  static int *partTable;
  static int updates=0;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  
  if(firsttime){
    tempList=malloc(sizeof(int)* pos.nAtoms*nd->nMax);
	 tempListnMax=nd->nMax;
	 tempHead=malloc(sizeof(int)* pos.nAtoms);
    partTable=malloc(sizeof(int)* pos.nAtoms);
	 firsttime=0;
  }

  if(tempListnMax<nd->nMax){	  /* if size of ngbrlist has been incewased outside this function.. */
    tempList=realloc(tempList,(size_t)sizeof(int)*pos.nAtoms*nd->nMax);
	 tempListnMax=nd->nMax;
	 printf("list increased to %d outside partNgbr\n",tempListnMax);
  }

  /*We are now updating the lists */
  updates++;
  /*init some values for checking of next update*/
  
  nd->forceUpdate=0;
  rcs2=nd->rcs*nd->rcs;  /* == (rc+rskin)*(rc+rskin);*/
 
  enforcePeriodicityPart(pos,par,partSkin,atomsInPartSkin);  /* make sure atoms are placed periodically */
  
  /* take copies of position so that we can calculate when we shoould update ngbrs the next time*/
  memset(partTable,0,sizeof(int)*pos.nAtoms);
  
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    partTable[i]=1; /*put in table which is in part */
    nd->prevx[i]=pos.x[i];
    nd->prevy[i]=pos.y[i];
    nd->prevz[i]=pos.z[i];
  }


  if(nd->volume!=par->volume){
    nd->volume=par->volume;
    cellInit(nd->cd,pos,par,0);
  }
  
  /*update cells first, this could still be part optimized..*/
  //cellUpdatePart(nd->cd,pos,par,partSkin,atomsInPartSkin);
  cellUpdate(nd->cd,pos,par);
  
  posInList=0;
  
  for(i=0;i<pos.nAtoms;i++){
    tempHead[i]=posInList;
    
    if(partTable[i]==0){ /*not in part, lets copy old results... */
      ndlistPos=nd->head[i]; /*where the ngbrs start*/
      while(nd->list[ndlistPos]!=-1) /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
		  tempList[posInList++]=nd->list[ndlistPos++];
      
    }
    else{
      n=nd->cd->cellList[i]; /*n is the cell of atom i*/
      for(k=0;k<27;k++) /*loop over ngbrcells*/
		  if(nd->cd->ngbrCells[n*27+k]>=0) /*if it is ok*/
			 for(j=nd->cd->head[nd->cd->ngbrCells[n*27+k]]; j >=0 ; j = nd->cd->list[j]){ /*loops over atoms in ngbrcell*/
				
				rx=pos.x[i]-pos.x[j];
				ry=pos.y[i]-pos.y[j];
				rz=pos.z[i]-pos.z[j];
				
				r2=LENGTH2(rx,ry,rz);
				
				if(r2>minhBox2){
				  PERIODIC(rx,ry,rz);
				  r2=LENGTH2(rx,ry,rz);
				}
				
				if( r2 <= rcs2 && j!=i){
				  tempList[posInList  ]=j;
				  posInList++;
				}
				
				/*if the ngbrs dont fit in list increase list size (-1 so that endmarker will always fit)*/
				if(posInList>=pos.nAtoms*nd->nMax-10){
				  nd->nMax=nd->nMax*1.5+1;
				  nd->list=realloc(nd->list,(size_t)sizeof(int)*pos.nAtoms*nd->nMax);
				  tempList=realloc(tempList,(size_t)sizeof(int)*pos.nAtoms*nd->nMax);
				  tempListnMax=nd->nMax;
				}
			 }
      /*leave marker for end of list*/
    } 
    tempList[posInList++]=-1;
	 
	 
  }
  
  /*copy result from temporary buffers to real buffers */
  memcpy(nd->list,tempList,(size_t)sizeof(int)*pos.nAtoms*nd->nMax);
  memcpy(nd->head,tempHead,(size_t)sizeof(int)*pos.nAtoms);
  
  
	 
  //  printf("nmax %d\n",nd->nMax);
  if(DEBUG && updates%1000==0){
    printf("ngbrPartlist updates %d\n",updates);
		 
  }
}
