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

#ifndef CELL_H
#define CELL_H
struct cellData
{
  int *list;     /*contains the lists of atoms in cells, this is a linked list*/
  int *head;     /* head[i] contains the index where the list of atoms that are in cell i
						  start in *list, this is in column major format..., that is
						  i=ix+nCellsX*(iy+nCellsY*iz);*/ 
  int *ngbrCells;/*contains a list of the ngbr cells of a certain cell*/

  int *cellList;     /*contains the cell in which each atom is*/


  int nCells;    /*total number of cells*/   
  int nCellsX;   /*number of cells in x          */
  int nCellsY;   /*                   y          */
  int nCellsZ;   /*               and z direction*/  
  double lCell;  /*the length of the sides of the cells */

};
void clearCells(struct cellData *cd); 
void cellInit(struct cellData *cd, struct systemPos pos,struct parameters *par,int firsttime);
void cellUpdate(struct cellData *cd, struct systemPos pos,struct parameters *par);
void copyCells(struct cellData *dest,struct cellData *source,int nAtoms);
void cellUpdatePart(struct cellData *cd, struct systemPos pos,struct parameters *par,int *part,int atomsInPart);
#endif
