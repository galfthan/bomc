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


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "shared.h"
#include "createLattice.h"
#include "initialize.h"
#include "bondList.h"
#include "randFunc.h"


/*The unctions in this module create different kinds of lattices */

void createZnBLattice(double boxx, double boxy, double boxz, int mat1, int mat2, 
		      double density, struct systemPos *pos,struct parameters *par)
{
  int i,j,k,l,iat;
  int ix,iy,iz;
  double xo,yo,zo,delta; /*xyzo is where the cell begins*/
  double a;
  double diabasis[8][3]={{0.0 , 0.0  , 0.0 },
			 {0.5 , 0.5  , 0.0 },
			 {0.5 , 0.0  , 0.5 },
			 {0.0 , 0.5  , 0.5 },
			 {0.25, 0.25 , 0.25}, 
			 {0.75, 0.75 , 0.25}, 
			 {0.75, 0.25 , 0.75}, 
			 {0.25, 0.75 , 0.75}}; 
  
  a=pow(8.0/density,1.0/3.0); /* this is the length of an unit cell */
  ix=(int)(boxx/a+0.5);			  /* make the closest possible amount of unit cells */
  iy=(int)(boxy/a+0.5);
  iz=(int)(boxz/a+0.5);
	
  initStructPos(pos,ix*iy*iz*8);  /* 8 atoms in a unit cell */
  

  par->box[0]=ix*a;
  par->box[1]=iy*a;
  par->box[2]=iz*a; 
  setSystemSize(par);

  delta=0.25*a;
  
  iat=-1;
  for (i=0;i<ix;i++) {
    xo=a*i+delta-par->hBox[0];
    for (j=0;j<iy;j++) {
      yo=a*j+delta-par->hBox[1];
      for (k=0;k<iz;k++) {
	zo=a*k+delta-par->hBox[2];
	for(l=0;l<8;l++){
	  iat++;
	  pos->x[iat]=xo+diabasis[l][0]*a;
	  pos->y[iat]=yo+diabasis[l][1]*a;
	  pos->z[iat]=zo+diabasis[l][2]*a;
	  if(l<4){
	    pos->aType[iat]=mat1;
	  }
	  else{
	    pos->aType[iat]=mat2;
	  }
			 
	}
			 
      }
    }
  }
    
  
}

void createIbcLattice(double boxx, double boxy, double boxz, int mat1, int mat2,double ratio,
		      double density, struct systemPos *pos,struct parameters *par)
{
  int i,j,k,l,iat;
  int ix,iy,iz;
  double xo,yo,zo,delta; /*xyzo is where the cell begins*/
  double a;
  double  m1Faces[4][3]={{0.0 , 0.0  , 0.0 },
			 {0.5 , 0.5  , 0.0 },
			 {0.5 , 0.0  , 0.5 },
			 {0.0 , 0.5  , 0.5 }};
  
  double m1Internal[4][3]={{0.25, 0.25 , 0.25}, 
			   {0.75, 0.75 , 0.25}, 
			   {0.75, 0.25 , 0.75}, 
			   {0.25, 0.75 , 0.75}};
  
  double m2basis[16][3];
  int m1Atoms,m2Atoms;
  int includeM2Atom[16];
  
  /*first put the m2 atoms for first internal m1*/
  for(i=0;i<4;i++)
    for (k=0;k<3;k++)  /*loop over xyz component*/
      m2basis[i][k]=0.5*(m1Internal[0][k]+m1Faces[i][k]);
  
  /*update the rest of o-basis based on first part*/
  for (i=1;i<4;i++) /*loop over internal si atoms*/
    for (j=0;j<4;j++) /*loop over oxygen ngbrs of si atom i*/
      for (k=0;k<3;k++)  /*loop over xyz component*/
	m2basis[i*4+j][k]=m2basis[j][k]+(m1Internal[i][k]-m1Internal[0][k]);
  
  /* these are the nymber of atoms per unitcell */
  m1Atoms=8;
  m2Atoms=(int)(m1Atoms*ratio);
  a=pow((m1Atoms+m2Atoms)/density,1.0/3.0); /* this is the length of an unit cell */
  
  ix=(int)(boxx/a+0.5);			  /* make the closest possible amount of unit cells */
  iy=(int)(boxy/a+0.5);
  iz=(int)(boxz/a+0.5);
  
  initStructPos(pos,ix*iy*iz*(m1Atoms+m2Atoms));
  par->box[0]=ix*a;
  par->box[1]=iy*a;
  par->box[2]=iz*a; 
  setSystemSize(par);

  /* here we create includeM2Atom[16];  which is a mask for the m2basis, using this we can include the correct number of m2 atoms in the system */
  for(i=0;i<16;i++) includeM2Atom[i]=0;
  i=0;								  
  while(i< m2Atoms){
    j=randNum(0)*16;
    if(includeM2Atom[j]==0){
      includeM2Atom[j]=1;
      i++;
    }
  }

  delta=0.25*a;
  iat=-1;
  for (i=0;i<ix;i++) {
    xo=a*i+delta-par->hBox[0];
    for (j=0;j<iy;j++) {
      yo=a*j+delta-par->hBox[1];
      for (k=0;k<iz;k++) {
	zo=a*k+delta-par->hBox[2];
	for(l=0;l<4;l++){
	  iat++;
	  pos->x[iat]=xo+m1Internal[l][0]*a;
	  pos->y[iat]=yo+m1Internal[l][1]*a;
	  pos->z[iat]=zo+m1Internal[l][2]*a;
	  pos->aType[iat]=mat1;
	  iat++;
	  pos->x[iat]=xo+m1Faces[l][0]*a;
	  pos->y[iat]=yo+m1Faces[l][1]*a;
	  pos->z[iat]=zo+m1Faces[l][2]*a;
	  pos->aType[iat]=mat1;
	}
	for(l=0;l<16;l++) if(includeM2Atom[l]){
	  iat++;
	  pos->x[iat]=xo+m2basis[l][0]*a;
	  pos->y[iat]=yo+m2basis[l][1]*a;
	  pos->z[iat]=zo+m2basis[l][2]*a;
	  pos->aType[iat]=mat2;
			 
	}

		  
      }
    }
  }
  
  
}

void createRandLattice(double boxx,double boxy,double boxz, double density, int mat1, int mat2, double proportion,
		       struct systemPos *pos,struct parameters *par)
{
  int i;
  int nAtomsMat1;
  int nAtomsMat2;
  

  par->box[0]=boxx;
  par->box[1]=boxy;
  par->box[2]=boxz; 
  setSystemSize(par);
  
  nAtomsMat1=proportion/(proportion+1)*density*par->volume;
  nAtomsMat2=1/(proportion+1)*density*par->volume; 
  initStructPos(pos,nAtomsMat1+nAtomsMat2); 
  
  for (i=0;i< nAtomsMat1+nAtomsMat2;i++) {
    pos->x[i]=randNum(0)*boxx-boxx/2.0;
    pos->y[i]=randNum(0)*boxy-boxy/2.0;
    pos->z[i]=randNum(0)*boxz-boxz/2.0;
    if(i<nAtomsMat1)
      pos->aType[i]=mat1;
    else
      pos->aType[i]=mat2;
  }
  
  
  
}
