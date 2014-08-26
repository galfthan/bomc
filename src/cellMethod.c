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
#include <limits.h>
#include <memory.h>

#include "shared.h"
#include "cellMethod.h"

/*this is a module that divides the atoms into cells. It can be used to speed up the calculation of
  ngbrlists (or used directly in potentials if one doesn't want to use ngbrlists). To use this module 
  one must first create a celldata structure and initialize it with cellInit(). After that one must
  update it each simulationstep before it is used with cellUpdate()
  
  struct cellData
  int *list;       contains the lists of atoms in cells, this is a linked list
  int *head;       head[i] contains the index where the list of atoms that are in cell i
                    start in *list, this is in column major format..., that is
						  i=ix+nCellsX*(iy+nCellsY*iz); 
  int *cellList;   contains the cell in which each atom is

  int *ngbrCells;  contains a list of the ngbr cells of a certain cell
  int nCells;      total number of cells   
  int nCellsX;     number of cells in x          
  int nCellsY;                        y          
  int nCellsZ;                    and z direction  
  double lCell;    the length of the sides of the cells 

  To get the atoms in cell i
  
  i=ix+nCellsX*(iy+nCellsY*iz)                     <- cell in pos ix,iy,iz
  for(atom=head[n]; atom >=0 ; atom = list[atom])  <- first atom in head[atom], after that get next with atom = list[atom]
                                                      -1 marks the end of list

*/
void clearCells(struct cellData *cd){
  free(cd->list);
  free(cd->head);
  free(cd->ngbrCells);
  free(cd->cellList);
   
}

void copyCells(struct cellData *dest,struct cellData *source,int nAtoms){
  if(dest->nCells!=source->nCells){
    dest->nCells=source->nCells;
    dest->head=realloc(dest->head,sizeof(int)*source->nCells);
    dest->ngbrCells=realloc(dest->ngbrCells,sizeof(int)*27*source->nCells);
  }
  
  memcpy(dest->list,source->list,sizeof(int)*nAtoms);
  memcpy(dest->head,source->head,sizeof(int)*source->nCells);
  memcpy(dest->ngbrCells,source->ngbrCells,sizeof(int)*27*source->nCells);
  memcpy(dest->cellList,source->cellList,sizeof(int)*nAtoms);
  dest->lCell=source->lCell;
  dest->nCellsX =source->nCellsX ;
  dest->nCellsY =source->nCellsY ;
  dest->nCellsZ =source->nCellsZ ;
}


void cellInit(struct cellData *cd, struct systemPos pos,struct parameters *par,int firsttime)
{
	int i;
	int nx,ny,nz,dx,dy,dz; 
	int ncx,ncy,ncz; 
	double lCelli;  /*inverse of cell size*/
	lCelli = 1.0/cd->lCell;

	/*we calculate the number of cells*/
	cd->nCellsX =(int)(par->box[0]*lCelli)  ; 
	cd->nCellsY =(int)(par->box[1]*lCelli)  ;
	cd->nCellsZ =(int)(par->box[2]*lCelli)  ;
	cd->nCells = cd->nCellsX*cd->nCellsY*cd->nCellsZ; 
  
	if(cd->nCells ==0){
	  printf("Too small system compared to skin+cutoff radius, trying to continue...box=%g %g %g \n" ,par->box[0],par->box[1],par->box[2]);
	  
	  cd->nCells=1;
	  cd->nCellsZ=1;
	  cd->nCellsY=1;
	  cd->nCellsX=1;
	}
	
	if(firsttime){
	  cd->list=NULL;
	  cd->head=NULL;
	  cd->ngbrCells=NULL;
	  cd->cellList=NULL;
	}	
	

	cd->list=realloc(cd->list,sizeof(int)*pos.nAtoms);
	cd->head=realloc(cd->head,sizeof(int)*cd->nCells);
	cd->ngbrCells=realloc(cd->ngbrCells,sizeof(int)*27*cd->nCells);
	cd->cellList=realloc(cd->cellList,sizeof(int)*pos.nAtoms);

	for(i=0;i<cd->nCells;i++)
	  cd->head[i]=-1;
	
	/*calculating ngbr cells, happens only ones, can be a bit slow..*/
	
	for(nx=0;nx<cd->nCellsX;nx++) /*loop over cells*/
	  for(ny=0;ny<cd->nCellsY;ny++)
      for(nz=0;nz<cd->nCellsZ;nz++){
		  
		  int count=0; /*the number of the ngbr cell we calculate 0...26*/ 
		  int iCell; /*index of cell*/
		  int niCell; /*index of ngbr cell*/
		  
		  iCell = nx+cd->nCellsX*ny+cd->nCellsX*cd->nCellsY*nz;
		  
		  for(dx=-1;dx<2;dx++)  /*loop over ngbrs of cells*/
			 for(dy=-1;dy<2;dy++)
				for(dz=-1;dz<2;dz++){
				  
				  /*take periodicity into account*/
				  /*ncx is nx+dx but with periodicities taken into account*/
				  
				  niCell=0;
				  if(par->periodic[0])
					 ncx=(nx+dx<0)?(cd->nCellsX-1):( (nx+dx>=cd->nCellsX)?(0):(nx+dx) );
				  else
					 ncx=nx+dx;
				  if(par->periodic[1])
					 ncy=(ny+dy<0)?(cd->nCellsY-1):( (ny+dy>=cd->nCellsY)?(0):(ny+dy) );
				  else
					 ncy=ny+dy;
				  
				  if(par->periodic[2])
					 ncz=(nz+dz<0)?(cd->nCellsZ-1):( (nz+dz>=cd->nCellsZ)?(0):(nz+dz) );
				  else
					 ncz=nz+dz;
				  
				  /*if not allowed*/
				  if(ncx<0 || ncx>=cd->nCellsX ||
					  ncy<0 || ncy>=cd->nCellsY ||
					  ncz<0 || ncz>=cd->nCellsZ)
					 niCell=-1; 
				  else 					 
					 niCell=ncx+cd->nCellsX*(ncy+cd->nCellsY*ncz);
				  
				  /*lets check if we have already put this cell in list
					 We have to make this so that atoms arent counted twice  
					 if we have less than 3 cells in some dimension*/
				  
				  for(i=0;i<count;i++)
					 if(niCell==cd->ngbrCells[iCell*27+i])
						niCell=-1;
				  
				  /*ok lets add to ngbr  list*/
				  cd->ngbrCells[iCell*27+count]=niCell;
				  count++;
				  
				}
      }
  
  
}

void cellUpdate(struct cellData *cd, struct systemPos pos,struct parameters *par){
  
  int i;
  int iCell,ixCell,iyCell,izCell;
  double lCelli;

  

  lCelli = 1.0/cd->lCell;
  
  /* sorting of the atoms into cells and making the linked list*/
  
  
  for(i=0;i<cd->nCells;i++)
	 cd->head[i]=-1;
  for(i=0;i<pos.nAtoms;i++){
	 ixCell=(pos.x[i] + par->hBox[0])*lCelli;
	 iyCell=(pos.y[i] + par->hBox[1])*lCelli;
	 izCell=(pos.z[i] + par->hBox[2])*lCelli;
	 
	 /*we cant have halffilled cells at the border=> cells that should be in halffilled cell umber cd->nCellsX  is really in the previous*/ 
	 ixCell=(ixCell==cd->nCellsX)?cd->nCellsX-1:ixCell;
	 iyCell=(iyCell==cd->nCellsY)?cd->nCellsY-1:iyCell;
	 izCell=(izCell==cd->nCellsZ)?cd->nCellsZ-1:izCell;
	 
	 
	 /*now check for errors */
		
	 /*sometimes if the system is not properly forced to being periodic the cellnumber might e negative =>not allowed,increse these with nCells*/
	
	 
	 ixCell=(ixCell<0)?0:ixCell;
	 iyCell=(iyCell<0)?0:iyCell;
	 izCell=(izCell<0)?0:izCell;
	 
	 ixCell =(ixCell>cd->nCellsX)?cd->nCellsX-1:ixCell;
	 iyCell =(iyCell>cd->nCellsY)?cd->nCellsY-1:iyCell;
	 izCell =(izCell>cd->nCellsZ)?cd->nCellsZ-1:izCell;
	 
	 iCell = ixCell+cd->nCellsX *(iyCell+cd->nCellsY*izCell);
	 
	 cd->list[i]=cd->head[iCell];
	 cd->head[iCell] = i; 
	 
	 cd->cellList[i]=iCell;
		
	 if(iCell<0 || iCell>= cd->nCells){
		printf("problem with cellmethod");
		exit(0);
	 }
	 
  }
	 
  
	 
}
