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


#include "cellMethod.h"
#ifndef NGBR_H
#define NGBR_H
struct ngbrData
{
  int *list;    /*contains the ngbrlists for the atoms, each separate list is ended with a -1 to mark the end*/
  int *head;    /*contains the starting position of the ngbrlist on list*/  
  double *prevx,*prevy,*prevz; /*contains the positions the atoms had last time  the lists where, used to check
											if we need to update the ngbrlist*/
  double volume;
  int nMax;            /*maximum number of ngbrs, it is automaically increased if neede*/
  double rcs;          /*rcut +skin,defines the maximum distance of the ngbrs from the atom*/
  struct cellData *cd; /*the data for the cellmethod which is used to speed up the ngbr calculation*/
  int forceUpdate;     /*if 1, the ngbrs are updated without checking if it is neccessary*/
};
#define INITLISTPOS(head,atomi,listpos) (listpos=head[atomi])
#define NEXTNGBR(list,ngbr,listpos) ((ngbr)=(list)[listpos++])
#define LOOPNGBR(list,ngbr,listpos) while(NEXTNGBR(list,ngbr,listpos)!=-1)



void clearNgbrs(struct ngbrData *nd);
void initNgbrs(struct ngbrData *nd, struct systemPos pos,struct parameters *par,double rcut);
void updateNgbrs(struct ngbrData *nd,struct systemPos pos,struct parameters *par);
void updateNgbrsPart(struct ngbrData *nd,struct systemPos pos,struct parameters *par,int *partSkin,int atomsInPartSkin);
int checkNgbrUpdatePart(struct ngbrData *nd, struct systemPos pos,struct parameters *par,int *partSkin,int atomsInPartSkin);
int checkNgbrUpdateAll(struct ngbrData *nd, struct systemPos pos,struct parameters *par);
void copyNgbrs(struct ngbrData *dest,struct ngbrData *source,int nAtoms);


#endif
