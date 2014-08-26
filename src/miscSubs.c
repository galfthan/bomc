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
#include <string.h>

#include "shared.h"
#include "parallel.h"
#include "randFunc.h"
#include "initialize.h"
#include "potentialWww.h"
#include "writeMat.h"
#include "fileio.h"
#include "conjgrad.h"

void calculateBulkModulus(struct systemPos pos,struct parameters *par){
  double oldBox[3];
  double scale;
  int i,j;
  double pot,potmin,vol;
  
  memcpy(pos.prevx,pos.x,sizeof(double)*3*pos.nAtoms);
  
  printf("scale volume scalepot minimizedpot\n");
  for(i=0;i<3;i++)
	 oldBox[i]=par->box[i];

  for(j=0;j<100;j++){
	 scale=1.0+(j-50)/500.0;
	 
	 for(i=0;i<3;i++)
		par->box[i]=oldBox[i]*scale;
    setSystemSize(par);
	 pot=wwwPot(par,pos,0)*TOEV;

    for(i=0;i<3*pos.nAtoms;i++)				  /* scale atom positions */
		pos.x[i]=pos.prevx[i]*scale;
	 vol=par->box[0]*par->box[1]*par->box[2];
	 potmin=fullCg(pos,par,0,0,1e10,wwwPot)*TOEV; /*optimization not allowed, accPotVal has a arbitrary value */
	 printf("%g %g %g %g\n",scale,vol,pot,potmin);	 
	 fflush(stdout);
  }
}

void checkForceCalc(struct systemPos pos,struct parameters *par){
  
  int steps=10000000;
  int i,j;
  double maxStep=0.0001;
  double orgPot,pot;
  double *drx,*dry,*drz;
  double integral=0.0;
  
  drx=malloc(sizeof(double)*pos.nAtoms);
  dry=malloc(sizeof(double)*pos.nAtoms);
  drz=malloc(sizeof(double)*pos.nAtoms);

  orgPot=wwwPot(par,pos,1);

  for(i=0;i<steps;i++){
	 j=randNum(1)*pos.nAtoms;
		drx[j]=j/pos.nAtoms*2.0*(0.5-randNum(1))*maxStep;
		dry[j]=2.0*(0.5-randNum(1))*maxStep;
		drz[j]=2.0*(0.5-randNum(1))*maxStep;
		pos.x[j]+=drx[j];
		pos.y[j]+=dry[j];
		pos.z[j]+=drz[j];

	  integral+=0.5*(drx[j]*pos.xa[j]+dry[j]*pos.ya[j]+drz[j]*pos.za[j]);
	
	  pot=wwwPot(par,pos,1);
	  integral+=0.5*(drx[j]*pos.xa[j]+dry[j]*pos.ya[j]+drz[j]*pos.za[j]);
	  if(i%50==0)
		 printf("%d    deltapot=%g force dot dx=%g  err=%g \n",i,(pot-orgPot)*TOEV/pos.nAtoms,-integral*TOEV/pos.nAtoms,((pot-orgPot)+integral)/(pot-orgPot+1e-30)); 
  }
  free(drx);
  free(drz);
  free(dry);
}

void enforcePeriodicity(struct systemPos pos,struct parameters *par){
  double max[3];
  double min[3];
  int c,i;


  /*take periodicity into account, if not periodic make sure system size is correct*/
  /*if the system is periodic we assume the partcles are between -halfbox...halfbox */
  
  for(i=0;i<pos.nAtoms;i++){
    if(par->periodic[0]){
      MAKEPERIODIC_DEF_DIM(pos.x[i],-par->hBox[0],par->hBox[0],par->box[0]);
    }
    else 
      if (pos.x[i] < min[0])
	min[0] = pos.x[i];
      else if (pos.x[i] >max[0])
	max[0] = pos.x[i];
    
    if(par->periodic[1]){
      MAKEPERIODIC_DEF_DIM(pos.y[i],-par->hBox[1],par->hBox[1],par->box[1]);
    }
    else 
      if (pos.y[i] < min[1])
	min[1] = pos.y[i];
      else if (pos.y[i] >max[1])
	max[1] = pos.y[i];
    
    if(par->periodic[2]){
      MAKEPERIODIC_DEF_DIM(pos.z[i],-par->hBox[2],par->hBox[2],par->box[2]);
    }	 
    else 
      if (pos.z[i] < min[2])
	min[2] = pos.z[i];
      else if (pos.z[i] >max[2])
	max[2] = pos.z[i];
  }
  
  
  /*if not periodic update the size of the sys*/ 
  
  /*
    for(c=0;c<3;c++)
    if(!par->periodic[c]){
    par->box[c]=2.0*MAX(fabs(max[c]),fabs(min[c]));
    setSystemSize(par);
    }
  */
}



void enforcePeriodicityPart(struct systemPos pos,struct parameters *par,int *partSkin, int atomsInPartSkin){
  int i,ips;


  /*take periodicity into account*/
  /*if the system is periodic we assume the partcles are between -halfbox...halfbox */
  for(ips=0;ips<atomsInPartSkin;ips++){
	 i=partSkin[ips];
    if(par->periodic[0])
      MAKEPERIODIC_DEF_DIM(pos.x[i],-par->hBox[0],par->hBox[0],par->box[0]);
    if(par->periodic[1])
      MAKEPERIODIC_DEF_DIM(pos.y[i],-par->hBox[1],par->hBox[1],par->box[1]);
	 if(par->periodic[2])
      MAKEPERIODIC_DEF_DIM(pos.z[i],-par->hBox[2],par->hBox[2],par->box[2]);
    }
  
}

void convertSeconds(int totSeconds,int *days,int *hours, int *minutes,int *seconds){
  	*days=totSeconds/86400;
	totSeconds-=86400*(*days);
	*hours=totSeconds/3600;
	totSeconds-=3600*(*hours);
	*minutes=totSeconds/60;
	totSeconds-=60*(*minutes);
	*seconds=totSeconds;
}


void calculateDynamicalMatrix(struct systemPos pos,struct parameters *par){
  double *dMatrix;
  double h=0.05;					  /* the testing displacement of the atoms */
  double **r,**a;                   /* the position matrix, a Nx3 matrix, built form the pos.x.. */
  double pot;
  int cgsteps;
   int i,j;
  int l,k;
  char filename[40];
  char matrixname[40];
  char temp[100];
  double back;
  /*
	 strcpy(filename,"vdos.par");
	 openParameterFile(filename);
	 getParValue("h",&h,"%lf"); 
	 getParValue("matrixname",&(matrixname[0]),"%s"); 
	 closeParameterFile();
  */
  printf("calc of dmatrix starting h=%g THIS FUNCTION HAS NOT BEEN TESTED CHECKED SINCE REWRITE!\nFirst CG part\n",h);
  
  dMatrix=malloc(sizeof(double)*POW2(3*pos.nAtoms));
  if(dMatrix==NULL){
	 printf("out of memory, tried to allocate %.2gMB\n",sizeof(double)*POW2(3*pos.nAtoms)/POW2(1024.0));
	 exit(0);
  }
  r=malloc(sizeof(double *)*3);
  a=malloc(sizeof(double *)*3);
  
  r[0]=pos.x;
  r[1]=pos.y;
  r[2]=pos.z;
  a[0]=pos.xa;
  a[1]=pos.ya;
  a[2]=pos.za;
  
  fullCg(pos,par,0,0,1e10,wwwPot); /*optimization not allowed, accPotVal has a arbitrary value */
  
  for(i=0;i<3*pos.nAtoms;i++)
	 for(j=0;j<3*pos.nAtoms;j++)
	   dMatrix[i*3*pos.nAtoms+j]=0.0;
  
  printf("calculating dynamical matrix, this is going to take a while...\n");
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    if(i>20 && i%(pos.nAtoms/20)==0)	printf("%.2g%%\n",100.0*(double)i/pos.nAtoms);
    for(k=0;k<3;k++){
	
      for(l=0;l<3;l++)
	for(j=0;j< pos.nAtoms;j++)
	  a[l][j]=0.0;

      back=r[k][i];				  
      r[k][i]=back+h;					  
      
      wwwPot(par,pos,1);	 
      for(j=0;j< pos.nAtoms;j++)
	for(l=0;l<3;l++){
	  dMatrix[(i*3+k)*3*pos.nAtoms+j*3+l]=a[l][j];
	  a[l][j]=0;
	}
		
      r[k][i]=back-h;
      wwwPot(par,pos,1);	 

      r[k][i]=back;					  /* now it is at original position */
		
      for(j=0;j< pos.nAtoms;j++)
	for(l=0;l<3;l++){
	  dMatrix[(i*3+k)*3*pos.nAtoms+j*3+l]=(a[l][j]-dMatrix[(i*3+k)*3*pos.nAtoms+j*3+l])/(2.0*h)*TOJOULE/(UNITDIST*UNITDIST);
	  a[l][j]=0;
	}
    }
  
  
  }
  
  writeMat("deriv_2",3*pos.nAtoms,3*pos.nAtoms,0,dMatrix,NULL);
  
  /* now write down the first derivatives to a matrix (should be zero i know.. */
  wwwPot(par,pos,1);	 
  for(i=0;i< pos.nAtoms;i++) /*loop over atoms i*/
	 for(k=0;k<3;k++)
		dMatrix[i*3+k]=a[k][i]*TOJOULE/(UNITDIST);
  writeMat("deriv_1",3*pos.nAtoms,1,0,dMatrix,NULL);

  /* now write down the positions of the atoms to a matrix */
  for(j=0;j< pos.nAtoms;j++){
	 dMatrix[j*3+0]=pos.x[j]*1e-10;
	 dMatrix[j*3+1]=pos.y[j]*1e-10;
	 dMatrix[j*3+2]=pos.z[j]*1e-10;
  }
  writeMat("pos",3*pos.nAtoms,3,0,dMatrix,NULL);

		
}

