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
#include <memory.h>
#include "shared.h"
#include "conjgrad.h"
#include "fileio.h"
#include "potentialVff.h"
#include "nrMinSubs.h"
#include "mcSubs.h"
#include "randFunc.h"
/* local global variables */

static double finalKt[3]; 
static double startKt[3]; 
static int finalCgStep;
static double scaleFact;
static int kMax;
static double wantedAccProb;


void initSimAnneal(struct systemPos pos){
  double temp;
  int error=0;
  
  openParameterFile("simAnneal.par");
  error+=getParValue("finalKt",&finalKt[AMORPH_FLAG],"%lf"); 
  error+=getParValue("startKt",&startKt[AMORPH_FLAG],"%lf"); 
  error+=getParValue("finalCgStep",&finalCgStep,"%d");
  error+=getParValue("scaleFact",&scaleFact,"%lf");
  error+=getParValue("steps_per_nAtoms_at_kT",&temp,"%lf");
  error+=getParValue("wanted_acceptance_probability",&wantedAccProb,"%lf");
  closeParameterFile();
  if(error){
	 printf("could not read needed parameters for simAnneal\n");
	 exit(0);
	 } 
	 startKt[AMORPH_FLAG]*=EV;
	 finalKt[AMORPH_FLAG]*=EV;
	
	 startKt[INTERFACE_FLAG]=startKt[AMORPH_FLAG];
	 finalKt[INTERFACE_FLAG]=finalKt[AMORPH_FLAG];
	 startKt[CRYST_FLAG]=0.1*startKt[AMORPH_FLAG];
	 finalKt[CRYST_FLAG]=0.1*finalKt[AMORPH_FLAG];
	 
	 kMax=(int)(temp*pos.nAtoms);

}


double simAnnealSphere(struct systemPos pos,struct parameters *par,double spherePos[3],double sphereR){
  int i,j,k;
  int atomi;
  double maxStep[3];
  double  accepted[3];
  double  tried[3];
  int flagVal;

  double newPos[3];
  double oldPos[3];
  double dr[3];
  double kT[3];

  double oldPot;
  double newPot;
  double tpot=0;
  int *sphereAtoms;
  int numOfSphereAtoms;
  double sphereR2;
  int localkMax;
  double r2,dx,dy,dz;

  sphereAtoms=malloc(sizeof(int)*pos.nAtoms);
  numOfSphereAtoms=0;
  sphereR2=POW2(sphereR);

  /* choose the atoms in the sphere */
  for(i=0;i<pos.nAtoms;i++){
	 dx=pos.x[i]-spherePos[0];
	 dy=pos.y[i]-spherePos[1];
	 dz=pos.z[i]-spherePos[2];
			  	  
	 r2=LENGTH2(dx,dy,dz);
	 if(r2>par->minhBox2){
		PERIODIC(dx,dy,dz);
		r2=LENGTH2(dx,dy,dz);
	 }
	 if(r2<sphereR2)
		sphereAtoms[numOfSphereAtoms++]=i;
  }

  localkMax=(int)((double)numOfSphereAtoms/pos.nAtoms*kMax); /* we make a local kmax that is scaled so that the same amount of 
																					 trial moves are made per atom as in the global simAnneal */
  resetOnePotVff(par,pos);		  /* makes sure that the ngrlists in potVff are correct */
  for(j=0;j<3;j++) kT[j]=startKt[j]/scaleFact; /* puts the starting cvalue of Kt */
  maxStep[AMORPH_FLAG]=0.5; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  maxStep[INTERFACE_FLAG]=0.5; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  maxStep[CRYST_FLAG]=0.1; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  
  for(j=0;j<3;j++) tried[j]=0;
  for(j=0;j<3;j++) accepted[j]=0;
  
  while(kT[0]>finalKt[0] || kT[1]>finalKt[1] || kT[2]>finalKt[2]){
	 	 
	 for(j=0;j<3;j++) kT[j]*=scaleFact; /* scale  temperature*/
	 
	 for(k=0;k<localkMax;k++){
		
		atomi=sphereAtoms[(int)(randNum(1)*numOfSphereAtoms)];
		flagVal=pos.aFlag[atomi];
		//while(flagVal==CRYST_FLAG && randNum(1)<1.0); /* bias toward the quartz,interface part */
	
		tried[flagVal]=1.0+0.995*tried[flagVal];
		accepted[flagVal]*=0.995;
		
		
		oldPos[0]=pos.x[atomi];
		oldPos[1]=pos.y[atomi];
		oldPos[2]=pos.z[atomi];
		 
		
		/* let`s calculate drx,dry,drz */
		oldPot=onePotVff(par,pos,atomi,oldPos);
		
		for(j=0;j<3;j++) dr[j]=maxStep[flagVal]*2.0*(0.5-randNum(1));
		
		/* we have calculated our dr, update newPos */
		for(j=0;j<3;j++) newPos[j]=oldPos[j]+dr[j];

		/* make sure its periodic */
		for(j=0;j<3;j++) if(par->periodic[j])
		  MAKEPERIODIC_DEF_DIM(newPos[j],-par->hBox[j],par->hBox[j],par->box[j]);
		
		/* new potential and force*/
		newPot=onePotVff(par,pos,atomi,newPos);
		 
			
		/* see if we allow it */
		if(newPot<oldPot || randNum(1)<exp(-(newPot-oldPot)/kT[flagVal])){ /*move accepted*/
		  pos.x[atomi]= newPos[0];
		  pos.y[atomi]= newPos[1];
		  pos.z[atomi]= newPos[2];
		  accepted[flagVal]+=1.0;
		}
		
		if(k%10 && accepted[flagVal]>tried[flagVal]*wantedAccProb)
		  maxStep[flagVal]*=1.0001;
		else
		  maxStep[flagVal]*=0.9999; 
		

	 }
	 
  }
  

  free(sphereAtoms);
  return forcesVff(par,pos,0);
}


double simAnneal(struct systemPos pos,struct parameters *par){
  int i,j,k;
  int atomi;
  double maxStep[3];
  double  accepted[3];
  double  tried[3];
  int flagVal;


  double newPos[3];
  double oldPos[3];
  double dr[3];
  double kT[3];

  double oldPot;
  double newPot;
  double tpot=0;

  

  resetOnePotVff(par,pos);		  /* makes sure that the ngrlists in potVff are correct */
  
  
  
  for(j=0;j<3;j++) kT[j]=startKt[j]/scaleFact;
  
  
  maxStep[AMORPH_FLAG]=0.5; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  maxStep[INTERFACE_FLAG]=0.5; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  maxStep[CRYST_FLAG]=0.1; //hardcoded but shouldn't matter very much, it will soon adapt anyway..*/
  
  for(j=0;j<3;j++) tried[j]=0;
  for(j=0;j<3;j++) accepted[j]=0;
  
  while(kT[0]>finalKt[0] || kT[1]>finalKt[1] || kT[2]>finalKt[2]){
	 
	 // printf("doris: %g\n", forcesVff(par,pos,0));
	 
	 for(j=0;j<3;j++) kT[j]*=scaleFact;
	 
	 //printf("%g maxStep %g\n",par->kT*TOEV,maxStep);
	 for(k=0;k<kMax;k++){

		atomi=(int)(randNum(1)*pos.nAtoms);
		flagVal=pos.aFlag[atomi];
		//while(flagVal==CRYST_FLAG && randNum(1)<1.0); /* bias toward the quartz,interface part */
		
		tried[flagVal]=1.0+0.995*tried[flagVal];
		accepted[flagVal]*=0.995;
		
		
		oldPos[0]=pos.x[atomi];
		oldPos[1]=pos.y[atomi];
		oldPos[2]=pos.z[atomi];
		 
		
		/* let`s calculate drx,dry,drz */
		oldPot=onePotVff(par,pos,atomi,oldPos);
		
		for(j=0;j<3;j++) dr[j]=maxStep[flagVal]*2.0*(0.5-randNum(1));
		
		/* we have calculated our dr, update newPos */
		for(j=0;j<3;j++) newPos[j]=oldPos[j]+dr[j];

		/* make sure its periodic */
		for(j=0;j<3;j++) if(par->periodic[j])
		  MAKEPERIODIC_DEF_DIM(newPos[j],-par->hBox[j],par->hBox[j],par->box[j]);
		
		/* new potential and force*/
		newPot=onePotVff(par,pos,atomi,newPos);
		/* see if we allow it */
		if(newPot<oldPot || randNum(1)<exp(-(newPot-oldPot)/kT[flagVal])){ /*move accepted*/
		  pos.x[atomi]= newPos[0];
		  pos.y[atomi]= newPos[1];
		  pos.z[atomi]= newPos[2];
		  accepted[flagVal]+=1.0;
		}
		/*
		  if(  k%300==0) {
		  printf(" %f %f %f  end\n",TOEV*kT[flagVal],forcesVff(par,pos,0)*TOEV/pos.nAtoms, maxStep[flagVal]);	
		  fflush(stdout);
		  }
		*/
		if(k%10 && accepted[flagVal]>tried[flagVal]*wantedAccProb)
		  maxStep[flagVal]*=1.0001;
		else
		  maxStep[flagVal]*=0.9999; 
		
		
	 }
	 
  } 
  if(finalCgStep){
	 conjGrad(pos,par,&i,&newPot,forcesVff);
	 return newPot;
  }
  
  return forcesVff(par,pos,0);
}
