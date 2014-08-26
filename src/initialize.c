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
#include <string.h>
#include "shared.h"
#include "initialize.h"


double gauss(double sigma, double mean,long int *seed);
double lcg(long int *seed);

void setSystemSize(struct parameters *par){ /* if par.box is set it updates the rest of the parameters */

  par->hBox[0]=0.5*par->box[0];
  par->hBox[1]=0.5*par->box[1];
  par->hBox[2]=0.5*par->box[2];
 
  par->volume= par->box[0]* par->box[1]* par->box[2];
  par->minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]); /* kanske man alltid borde rakna det har... */
  par->minhBox2*=par->minhBox2; /*make sure its the ^2 .*/ 

}

void initStructParameters(struct parameters *par){
  par->periodic[0]=0;
  par->periodic[1]=0;
  par->periodic[2]=0;
  par->box[0]=-1;
  par->box[1]=-1;
  par->box[2]=-1;

  par->hBox[0]=-1;
  par->hBox[1]=-1;
  par->hBox[2]=-1;
  par->volume=-1;
  par->skin=0.2;
  
  /*atypes are:
	 0  1  2  3  4  5  6  7  8  9  10 11 12
	 Cu Si Ge C  Ga As In P  Al O  Ar Na Cl
  */ 

  strcpy(par->aName[0],"O "); /*important, the length of the name have to be 2, */
  strcpy(par->aName[1],"Si");
  par->aMass[0]=15.9994; /*in u*/
  par->aMass[1]=28.0855; /*in u*/

  par->kT=0.05*EV;

  par->interfaceSys=0;

  par->seed=malloc(sizeof(long));
  *(par->seed)=1;
}

void initStructPos(struct systemPos *pos ,int nAtoms){
  int i;
  pos->nAtoms=nAtoms;  /* 4 atoms in a unit cell in FCC */
  
  
  pos->x=malloc(sizeof(double)*pos->nAtoms*3);
  pos->y=pos->x+pos->nAtoms;
  pos->z=pos->x+2*pos->nAtoms;
  pos->prevx=malloc(sizeof(double)*pos->nAtoms*3);
  pos->prevy=pos->prevx+pos->nAtoms;
  pos->prevz=pos->prevx+2*pos->nAtoms;
  pos->xa=malloc(sizeof(double)*pos->nAtoms*3);
  pos->ya=pos->xa+pos->nAtoms;
  pos->za=pos->xa+2*pos->nAtoms;
  pos->enPot= malloc(sizeof(double)*pos->nAtoms);
  pos->prevEnPot= malloc(sizeof(double)*pos->nAtoms);
  pos->aType=malloc(sizeof(int)*pos->nAtoms);
  pos->aFlag=malloc(sizeof(long)*pos->nAtoms);
  

  if(pos->x==NULL  || pos->y==NULL  || pos->z==NULL || 
	  pos->xa==NULL || pos->ya==NULL || pos->za==NULL || 
	  pos->enPot==NULL || pos->aType==NULL){
	 printf("ERROR! could not allocate memory in initStructPos()\n");
	 exit(0);
  }

 for(i=0;i<pos->nAtoms;i++){
   
	 pos->x[i]=0.0;
	 pos->y[i]=0.0;
	 pos->z[i]=0.0;
	 pos->xa[i]=0.0;
	 pos->ya[i]=0.0;
	 pos->za[i]=0.0;
	 pos->enPot[i]=0.0;
	 pos->aType[i]=-1;
	 pos->aFlag[i]=0;
 }
  
  
}



double lcg(long int *seed){
  static long int a=69069,c=1,m=2147483647;
  static double rm=2147483647.0;
  
  *seed=(*seed*a+c)%m;
  return (double)*seed/rm;

}



double gauss(double sigma, double mean,long int *seed){
  
  double rgauss;
  double v1, v2;
  double r=2.0;
  while(r >= 1.0){
    v1=2.0*(double)(lcg(seed)-1.0);
    v2=2.0*(double)(lcg(seed)-1.0);
    r=v1*v1+v2*v2;
  }
  rgauss = v1*sqrt(-2.0*log(r)/r);
  rgauss = mean+rgauss*sigma;
  return rgauss; 
}



  

