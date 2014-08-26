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


#include "ngbrs.h"

struct bondList{
  int nAtoms;
  int nMax; /*max amount of ngbrs*/
  int *list;    /*contains the ngbrlists for the atoms, each separate list is ended with a -1 to mark the end*/
  int *head;    /*contains the starting position of the ngbrlist on list*/
  
  int *danglBonds; /*contains the amount of dangling bonds per atom*/
  int arraySize;

};

void deleteLongBonds(struct systemPos pos,struct parameters *par,double rmax);

void initBondList(struct systemPos pos);
void connectAtoms(struct systemPos pos,struct parameters *par,double rcut);

void getBondInfo(int **head,int **list);
void getBondList(struct bondList  **lblist);

void nonNnNgbrsPart(struct ngbrData *nd,int nAtoms,int *partSkin,int atomsInPartSkin);
void nonNnNgbrs(struct ngbrData *nd,int nAtoms);

int removeBond(int ai, int aj);
int createBond(struct systemPos *pos,int ai, int aj);

int areBonded(int ai,int aj);
int areBonded2(int ai,int aj);
int numBonds(int ai);

int getRandAtomWb(int minBonds);
int getRandAtomWdb(int minDanglBonds);
int getRandBondAtAtom(int ai);
int getNextNode(int a1, int a2);

