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
#include <limits.h>
#include <time.h>
#include "shared.h"
#include "fileio.h"
#include "ngbrs.h"
#include "potentialVff.h"
#include "mcSubs.h"
#include "miscSubs.h"
#include "bondList.h"
#include "parallel.h"

static int siType=1,oType=0;

struct potentialVff
{
  double repR0;
  double repR02;
  double repK;
  double rcut,rcut2;
  
  double r0[2][2];
  double cos0[2];
  double kStretch[2][2];
  double kBend[2][2][2];
  
  double subOxidePenalty[5];
   
  int includeRepPot;
  int includeVffPot;
  int includeSuboxidePot;
  int includeDanglbondPot;
  
};




#define KBEND(i,j,k) (vff.kBend[pos.aType[i]][pos.aType[j]][pos.aType[k]])
#define KSTRETCH(i,j) (vff.kStretch[pos.aType[i]][pos.aType[j]])
#define COS0(i) (vff.cos0[pos.aType[i]])
#define R0(i,j) (vff.r0[pos.aType[i]][pos.aType[j]])


static struct potentialVff vff;
static struct ngbrData nd;
static struct ngbrData ndLastStep;




void initPotentialVff(struct parameters *par,struct systemPos pos){
  int error,i,j,k;
  char filename[20];
  strcpy(filename,"potVff.par");
  openParameterFile(filename);
  error=0;
  
  error+=getParValue("includeVffPot",&vff.includeVffPot,"%d");
  error+=getParValue("includeRepPot",&vff.includeRepPot,"%d");
  error+=getParValue("includeSuboxidePot",&vff.includeSuboxidePot,"%d");
  error+=getParValue("includeDanglbondPot",&vff.includeDanglbondPot,"%d");


  error+=getParValue("kStretch_SiSi",&vff.kStretch[1][1],"%lf");
  error+=getParValue("kStretch_SiO",&vff.kStretch[1][0],"%lf");
 
  error+=getParValue("kBend_SiSiSi",&vff.kBend[1][1][1],"%lf");
  error+=getParValue("kBend_SiOSi",&vff.kBend[1][0][1],"%lf");
  error+=getParValue("kBend_OSiO",&vff.kBend[0][1][0],"%lf");
  
  
  error+=getParValue("r0_SiSi",&vff.r0[1][1],"%lf");
  error+=getParValue("r0_SiO",&vff.r0[1][0],"%lf");
  
  error+=getParValue("cos0_Si",&vff.cos0[1],"%lf");
  error+=getParValue("cos0_O",&vff.cos0[0],"%lf");
  
  
  error+=getParValue("repulsive_r0",&vff.repR0,"%lf");
  error+=getParValue("repulsive_k",&vff.repK,"%lf");

  error+=getParValue("SubOxidePenalty+0",&vff.subOxidePenalty[0],"%lf");
  error+=getParValue("SubOxidePenalty+1",&vff.subOxidePenalty[1],"%lf");
  error+=getParValue("SubOxidePenalty+2",&vff.subOxidePenalty[2],"%lf");
  error+=getParValue("SubOxidePenalty+3",&vff.subOxidePenalty[3],"%lf");
  error+=getParValue("SubOxidePenalty+4",&vff.subOxidePenalty[4],"%lf");

  closeParameterFile();
    
  vff.kBend[1][1][0]=sqrt(vff.kBend[0][1][0]*vff.kBend[1][0][1]);
  vff.kBend[0][1][1]=vff.kBend[1][1][0];
  
  
  /*make the matrices symmetric*/
  vff.kStretch[0][1]=vff.kStretch[1][0];
  vff.r0[0][1]=vff.r0[1][0];
  
  printf("potVff: include rep %d\n",vff.includeRepPot);

  
  
  /*change the units*/
  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      vff.kStretch[i][j]*=EV;
      for(k=0;k<2;k++)
	vff.kBend[i][j][k]*=EV;
    }
  for(i=0;i<5;i++)
    vff.subOxidePenalty[i]*=EV;
  
  vff.repK*=EV;
  vff.rcut=vff.repR0;
  vff.rcut2=POW2(vff.rcut);
  vff.repR02=POW2(vff.repR0);

 
  
  if(error!=0){
    printf("could not read the required paramters from file %s\n",filename);
    exit(0);
  }
    
  /*init ngbrs */
  initNgbrs(&nd,pos,par,vff.rcut);
  initNgbrs(&ndLastStep,pos,par,vff.rcut);
  
  updateNgbrs(&nd,pos,par);
  nonNnNgbrs(&nd,pos.nAtoms); 
  copyNgbrs(&ndLastStep,&nd,pos.nAtoms);
  
  
  if(getRank()==0) {
    printf("potential initialized\n");
	 
    printf("K_r  Si-Si: %g ev/A^2\n",vff.kStretch[1][1]*TOEV);
    printf("     O-O:   %g ev/A^2\n",vff.kStretch[0][0]*TOEV);
    printf("     Si-O:  %g ev/A^2\n",vff.kStretch[1][0]*TOEV);
    printf("     O-Si:  %g ev/A^2\n",vff.kStretch[0][1]*TOEV);
	 
    printf("r0   Si-Si: %g A\n",vff.r0[1][1]);
    printf("     O-O:   %g A\n",vff.r0[0][0]);
    printf("     Si-O:  %g A\n",vff.r0[1][0]);
    printf("     O-Si:  %g A\n\n",vff.r0[0][1]);
	 
    printf("theta0 Si:       %g \n",acos(vff.cos0[1])*180.0/M_PI);
    printf("       O:        %g \n",acos(vff.cos0[0])*180.0/M_PI);
	 
    printf("K_bend Si-Si-Si: %g ev\n",vff.kBend[1][1][1]*TOEV);
    printf("       O-Si-O:   %g ev\n",vff.kBend[0][1][0]*TOEV);
    printf("       Si-O-Si:  %g ev\n",vff.kBend[1][0][1]*TOEV);
    printf("       Si-Si-O:  %g ev\n",vff.kBend[1][1][0]*TOEV);
    printf("       O-Si-Si:  %g ev\n",vff.kBend[0][1][1]*TOEV);
    printf("       Si-O-O:   %g ev\n",vff.kBend[1][0][0]*TOEV);
    printf("       O-O-Si:   %g ev\n",vff.kBend[0][0][1]*TOEV);
    printf("       O-O-O:    %g ev\n",vff.kBend[0][0][0]*TOEV);

  }

}


double forcesVff(struct parameters *par,struct systemPos pos,     int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  int i,j,k;
  double rij,rij2;
  double rik,rik2;
  
  int numberOfOngbrs;
  
  int ngbrListPosj;
  int bondListPosj,bondListPosk;
  double forcej,forcek;
  double varTemp1,varTemp2,varTemp3,varTemp4;
  double cosijk,cosTerm;
  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;
  int bondij,bondik;
  static double *bondR;
  static double *bondR2;
  static double *bondDx;
  static double *bondDy;
  static double *bondDz;

  static int firsttime=1;
  struct bondList *blist;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  getBondList(&blist);


 
  
  
  /*update ngbrs*/
  if(vff.includeRepPot && checkNgbrUpdateAll(&nd,pos,par)){
    updateNgbrs(&nd,pos,par);
    nonNnNgbrs(&nd,pos.nAtoms); /*tkaes away first and scond ngbrs that are bonded*/
  }

  if(calcForces){ /*this is not good if we want to have more than one pot at the same time!!!*/
    memset(pos.xa,0,pos.nAtoms*sizeof(double));
    memset(pos.ya,0,pos.nAtoms*sizeof(double));
    memset(pos.za,0,pos.nAtoms*sizeof(double));
  }
  if(firsttime){
    firsttime=0;
    bondR=malloc(4*pos.nAtoms*sizeof(double));
    bondR2=malloc(4*pos.nAtoms*sizeof(double));
    bondDx=malloc(4*pos.nAtoms*sizeof(double));
    bondDy=malloc(4*pos.nAtoms*sizeof(double));
    bondDz=malloc(4*pos.nAtoms*sizeof(double));
  }
   
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    bondListPosj=blist->head[i]; /*where the ngbrs start*/
    while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=blist->list[bondListPosj]; /*j is the ngbr*/
      bondij=4*i+bondListPosj-blist->head[i];
      bondListPosj++;
      bondDx[bondij]=pos.x[j]-pos.x[i]; 
      bondDy[bondij]=pos.y[j]-pos.y[i]; 
      bondDz[bondij]=pos.z[j]-pos.z[i]; 
	  
      bondR2[bondij]=LENGTH2(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
      if(bondR2[bondij]>minhBox2){
	PERIODIC(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
	bondR2[bondij]=LENGTH2(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
      }
	  
      bondR[bondij]=sqrt(bondR2[bondij]);
    }
  } 
  totPot=0.0;
  /*calculate forces*/
  if(vff.includeVffPot)  /*dont know why one would want to exclude this part but...*/
    for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
      bondListPosj=blist->head[i]; /*where the ngbrs start*/
      numberOfOngbrs=0; 
      while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
	j=blist->list[bondListPosj]; /*j is the ngbr*/
    	bondij=4*i+bondListPosj-blist->head[i];
	bondListPosj++;
		
	/*	 calculate number of oxygen ngbrs */ 
	if(pos.aType[j]==oType)
	  numberOfOngbrs++;		  
		
		
	drx_ij=bondDx[bondij];
	dry_ij=bondDy[bondij];
	drz_ij=bondDz[bondij];
	rij2=bondR2[bondij];
	rij=bondR[bondij];
		
	/*first the pair part*/
	/*here we sum over  i=0...N-1, j>i*/
      
	if(j>i){ /*dont doublecount!! in pair part*/
	  poti=KSTRETCH(i,j)*POW2(R0(i,j)-rij);
	  totPot+=poti;
	  pos.enPot[i]=0.5*poti;
	  pos.enPot[j]=0.5*poti;

	  if(calcForces){ /*j=n*/
	    forcej=2.0*KSTRETCH(i,j)*(R0(i,j)-rij)/rij;
	    /*xyza is acceleration, NOT force*/
	    pos.xa[i] -= forcej*drx_ij;
	    pos.ya[i] -= forcej*dry_ij;
	    pos.za[i] -= forcej*drz_ij;
			 
	    pos.xa[j] += forcej*drx_ij;
	    pos.ya[j] += forcej*dry_ij;
	    pos.za[j] += forcej*drz_ij;
	  }
	}
		
      
	/*then the three body part */
	/*here we sum over i=0...N-1, j!=i, k!=i && k>j */
      
	bondListPosk=bondListPosj; /*bondListPosj was increased already erlier so we dont need +1*/
	while(blist->list[bondListPosk]!=-1){ /*get second ngbr k>j*/
	  k=blist->list[bondListPosk]; /*k is the ngbr*/
	  bondik=4*i+bondListPosk-blist->head[i];
	  bondListPosk++;
		  
	  drx_ik=bondDx[bondik];
	  dry_ik=bondDy[bondik];
	  drz_ik=bondDz[bondik];
	  rik2=bondR2[bondik];
	  rik=bondR[bondik]; 
		
	  cosijk=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/(rik*rij);
	  cosTerm=(cosijk-COS0(i));
	  poti=KBEND(j,i,k)*POW2(cosTerm);
	  pos.enPot[i]=poti;
	  totPot+=poti;

		  
	  if(calcForces){
	    varTemp1=-2.0*KBEND(j,i,k)*cosTerm;
	    varTemp2=varTemp1/(rik*rij);
	    varTemp3=varTemp1*cosijk/POW2(rij);
	    varTemp4=varTemp1*cosijk/POW2(rik);
			

	    forcej=varTemp2*drx_ik-varTemp3*drx_ij;
	    forcek=varTemp2*drx_ij-varTemp4*drx_ik;
	    pos.xa[j] += forcej;
	    pos.xa[k] += forcek;
	    pos.xa[i] += -forcej-forcek;
			 
	    forcej=varTemp2*dry_ik-varTemp3*dry_ij;
	    forcek=varTemp2*dry_ij-varTemp4*dry_ik;
	    pos.ya[j] += forcej;
	    pos.ya[k] += forcek;
	    pos.ya[i] += -forcej-forcek;
			 
	    forcej=varTemp2*drz_ik-varTemp3*drz_ij;
	    forcek=varTemp2*drz_ij-varTemp4*drz_ik;
	    pos.za[j] += forcej;
	    pos.za[k] += forcek;
	    pos.za[i] += -forcej-forcek;
	
	  }
		  
	}
      
      
      }
      if(vff.includeSuboxidePot)
	if(pos.aType[i]==siType){
	  poti=vff.subOxidePenalty[numberOfOngbrs];
	  pos.enPot[i]=poti;
	  totPot+=poti;
	}
	 
	 

    }
 

  /* dangl Bond energy */
  if(vff.includeDanglbondPot)
    for(i=0;i < pos.nAtoms;i++){ 
      if(pos.aType[i]==siType)
	poti=blist->danglBonds[i]*1.0*EV;
      else if(pos.aType[i]==oType)
	poti=blist->danglBonds[i]*4.0*EV;
	  
      pos.enPot[i]=poti;
      totPot+=poti;
    }
 
 
  /*calculate the pot part that makes sure that the atoms do not overlap if not bonded*/
  k=0;
 
  for(i=0;i < pos.nAtoms && vff.includeRepPot;i++){ /*loop over atoms i*/
    //for(i=pos.nAtoms-1;i>=0;i--){ /*loop over atoms i, from biggest to smallest*/
    ngbrListPosj=nd.head[i]; /*where the ngbrs start*/
    while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=nd.list[ngbrListPosj]; /*j is the ngbr*/
      ngbrListPosj++;
      if(i>j){ /*dont doublecount!! in pair part*/
		 
	drx_ij=pos.x[j]-pos.x[i]; 
	dry_ij=pos.y[j]-pos.y[i]; 
	drz_ij=pos.z[j]-pos.z[i];	
	rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		  
	if(rij2>minhBox2){
	  PERIODIC(drx_ij,dry_ij,drz_ij);
	  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
	}
	if( rij2 <vff.repR02 ){
	  //	 printf("atoms %d %d rep %g\n",i,j,sqrt(rij2));
			 
	  k++;			 
	  rij=sqrt(rij2);
	  poti=vff.repK*POW2(vff.repR0-rij);
	  //poti=vff.repK*(vff.repR02-rij2);
	  totPot+=poti;
	  pos.enPot[i]+=0.5*poti;
	  pos.enPot[j]+=0.5*poti;
	  if(calcForces){ /*j=n*/
	    forcej=2.0*vff.repK*(vff.repR0-rij)/rij;
	    //forcej=2.0*vff.repK;
	    /*xyza is acceleration, NOT force*/
	    pos.xa[i] -= forcej*drx_ij;
	    pos.ya[i] -= forcej*dry_ij;
	    pos.za[i] -= forcej*drz_ij;
				
	    pos.xa[j] += forcej*drx_ij;
	    pos.ya[j] += forcej*dry_ij;
	    pos.za[j] += forcej*drz_ij;
	  }	
	}
      }
    }
  }


  
  firsttime=0;
  return totPot;
}



double partRepulsive(struct parameters *par,struct systemPos pos,int *partSkin, int atomsInPartSkin,int calcForces){
  int i,ips,j;
  int ngbrListPosj;
  double poti,totPot=0.0;
  double forcej;
  double drx_ij,dry_ij,drz_ij;
  double rij,rij2;
	
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  
  /*update ngbrs*/
  if(checkNgbrUpdatePart(&nd,pos,par,partSkin,atomsInPartSkin)){
    //	 updateNgbrs(&nd,pos,par);
    updateNgbrsPart(&nd,pos,par,partSkin,atomsInPartSkin);
    nonNnNgbrsPart(&nd,pos.nAtoms,partSkin,atomsInPartSkin); /*tkaes away first and scond ngbrs that are bonded*/
  }
  
  
  /*calculate the pot part that makes sure that the atoms do not overlap if not bonded*/
  
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
	 
    //for(i=pos.nAtoms-1;i>=0;i--){ /*loop over atoms i, from biggest to smallest*/
    ngbrListPosj=nd.head[i]; /*where the ngbrs start*/
    while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=nd.list[ngbrListPosj]; /*j is the ngbr*/
      ngbrListPosj++;
      if(i>j){ /*dont doublecount!! in pair part*/
		  
	drx_ij=pos.x[j]-pos.x[i]; 
	dry_ij=pos.y[j]-pos.y[i]; 
	drz_ij=pos.z[j]-pos.z[i];	
	rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		  
	if(rij2>minhBox2){
	  PERIODIC(drx_ij,dry_ij,drz_ij);
	  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
	}
	if( rij2 <vff.repR02 ){
	  //	 printf("atoms %d %d rep %g\n",i,j,sqrt(rij2));
			 
			 
	  rij=sqrt(rij2);
	  poti=vff.repK*POW2(vff.repR0-rij);
	  //poti=vff.repK*(vff.repR02-rij2);
	  totPot+=poti;
	  pos.enPot[i]+=0.5*poti;
	  pos.enPot[j]+=0.5*poti;
	  if(calcForces){ /*j=n*/
	    forcej=2.0*vff.repK*(vff.repR0-rij)/rij;
	    //forcej=2.0*vff.repK;
	    /*xyza is acceleration, NOT force*/
	    pos.xa[i] -= forcej*drx_ij;
	    pos.ya[i] -= forcej*dry_ij;
	    pos.za[i] -= forcej*drz_ij;
				
	    pos.xa[j] += forcej*drx_ij;
	    pos.ya[j] += forcej*dry_ij;
	    pos.za[j] += forcej*drz_ij;
	  }	
	}
      }
    }
  }
  return totPot;
}


double partSuboxide(struct parameters *par,struct systemPos pos,int *part, int atomsInPart){
  int i,ip,j;
  int numberOfOngbrs,bondListPosj;
  double poti,totPot=0.0;
  struct bondList *blist;
  getBondList(&blist);
  
  for(ip=0;ip<atomsInPart;ip++){
    i=part[ip];
    if(pos.aType[i]==siType){
      bondListPosj=blist->head[i]; /*where the ngbrs start*/
      numberOfOngbrs=0; 
      while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
	j=blist->list[bondListPosj]; /*j is the ngbr*/
	bondListPosj++;
	if(pos.aType[j]==oType)
	  numberOfOngbrs++;		  
      }
		
		
      poti=vff.subOxidePenalty[numberOfOngbrs];
      pos.enPot[i]=poti;
      totPot+=poti;
    } 
  }
  return totPot;
}

double partDanglBonds(struct parameters *par,struct systemPos pos,int *part, int atomsInPart){
  int i,ip;
  double poti,totPot=0.0;
  struct bondList *blist;
  getBondList(&blist);
  
  for(ip=0;ip<atomsInPart;ip++){
    i=part[ip];
    if(pos.aType[i]==siType)
      poti=blist->danglBonds[i]*1.0*EV;
    else if(pos.aType[i]==oType)
      poti=blist->danglBonds[i]*4.0*EV;
	 
    pos.enPot[i]=poti;
    totPot+=poti;
  }
  return totPot;
}
double partKeating(struct parameters *par,struct systemPos pos,int *partSkin,int atomsInPartSkin,int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  int i,ips,j,k;
  double rij,rij2;
  double rik,rik2;
  
 
  int bondListPosj,bondListPosk;
  double forcej,forcek;
  double varTemp1,varTemp2,varTemp3,varTemp4;
  double cosijk,cosTerm;
  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;
  int bondij,bondik;
  static double *bondR;
  static double *bondR2;
  static double *bondDx;
  static double *bondDy;
  static double *bondDz;

  static int firsttime=1;
  struct bondList *blist;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  getBondList(&blist);
  

 
    
  if(firsttime){
    firsttime=0;
    bondR=malloc(4*pos.nAtoms*sizeof(double));
    bondR2=malloc(4*pos.nAtoms*sizeof(double));
    bondDx=malloc(4*pos.nAtoms*sizeof(double));
    bondDy=malloc(4*pos.nAtoms*sizeof(double));
    bondDz=malloc(4*pos.nAtoms*sizeof(double));
  }
  
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    bondListPosj=blist->head[i]; /*where the ngbrs start*/
    while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=blist->list[bondListPosj]; /*j is the ngbr*/
      bondij=4*ips+bondListPosj-blist->head[i];
      bondListPosj++;
      bondDx[bondij]=pos.x[j]-pos.x[i]; 
      bondDy[bondij]=pos.y[j]-pos.y[i]; 
      bondDz[bondij]=pos.z[j]-pos.z[i]; 
		
      bondR2[bondij]=LENGTH2(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
      if(bondR2[bondij]>minhBox2){
	PERIODIC(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
	bondR2[bondij]=LENGTH2(bondDx[bondij],bondDy[bondij],bondDz[bondij]);
      }
		
      bondR[bondij]=sqrt(bondR2[bondij]);
    }
  }
 
  totPot=0.0;
  /*calculate forces*/
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    bondListPosj=blist->head[i]; /*where the ngbrs start*/
    while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=blist->list[bondListPosj]; /*j is the ngbr*/
      bondij=4*ips+bondListPosj-blist->head[i];
      bondListPosj++;
		
		
      drx_ij=bondDx[bondij];
      dry_ij=bondDy[bondij];
      drz_ij=bondDz[bondij];
      rij2=bondR2[bondij];
      rij=bondR[bondij];
		
      /*first the pair part*/
      /*here we sum over  i=0...N-1, j>i*/
      
      if(j>i){ /*dont doublecount!! in pair part*/
	poti=KSTRETCH(i,j)*POW2(R0(i,j)-rij);
	totPot+=poti;
	pos.enPot[i]=0.5*poti;
	pos.enPot[j]=0.5*poti;
			
	if(calcForces){ /*j=n*/
	  forcej=2.0*KSTRETCH(i,j)*(R0(i,j)-rij)/rij;
	  /*xyza is acceleration, NOT force*/
	  pos.xa[i] -= forcej*drx_ij;
	  pos.ya[i] -= forcej*dry_ij;
	  pos.za[i] -= forcej*drz_ij;
			  
	  pos.xa[j] += forcej*drx_ij;
	  pos.ya[j] += forcej*dry_ij;
	  pos.za[j] += forcej*drz_ij;
	}
      }
		 
		 
      /*then the three body part */
      /*here we sum over i=0...N-1, j!=i, k!=i && k>j */
		 
      bondListPosk=bondListPosj; /*bondListPosj was increased already erlier so we dont need +1*/
      while(blist->list[bondListPosk]!=-1){ /*get second ngbr k>j*/
	k=blist->list[bondListPosk]; /*k is the ngbr*/
	bondik=4*ips+bondListPosk-blist->head[i];
	bondListPosk++;
		  
	drx_ik=bondDx[bondik];
	dry_ik=bondDy[bondik];
	drz_ik=bondDz[bondik];
	rik2=bondR2[bondik];
	rik=bondR[bondik]; 
		
	cosijk=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/(rik*rij);
	cosTerm=(cosijk-COS0(i));
	poti=KBEND(j,i,k)*POW2(cosTerm);
	pos.enPot[i]=poti;
	totPot+=poti;

		  
	if(calcForces){
	  varTemp1=-2.0*KBEND(j,i,k)*cosTerm;
	  varTemp2=varTemp1/(rik*rij);
	  varTemp3=varTemp1*cosijk/rij2;
	  varTemp4=varTemp1*cosijk/rik2;
			

	  forcej=varTemp2*drx_ik-varTemp3*drx_ij;
	  forcek=varTemp2*drx_ij-varTemp4*drx_ik;
	  pos.xa[j] += forcej;
	  pos.xa[k] += forcek;
	  pos.xa[i] += -forcej-forcek;
			 
	  forcej=varTemp2*dry_ik-varTemp3*dry_ij;
	  forcek=varTemp2*dry_ij-varTemp4*dry_ik;
	  pos.ya[j] += forcej;
	  pos.ya[k] += forcek;
	  pos.ya[i] += -forcej-forcek;
			 
	  forcej=varTemp2*drz_ik-varTemp3*drz_ij;
	  forcek=varTemp2*drz_ij-varTemp4*drz_ik;
	  pos.za[j] += forcej;
	  pos.za[k] += forcek;
	  pos.za[i] += -forcej-forcek;
	
	}
		  
      }
      
      
    }
  }
  return totPot;
}
double partKeatingNotable(struct parameters *par,struct systemPos pos,int *partSkin,int atomsInPartSkin,int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  int i,ips,j,k;
  double rij,rij2;
  double rik,rik2;
  
 
  int bondListPosj,bondListPosk;
  double forcej,forcek;
  double varTemp1,varTemp2,varTemp3,varTemp4;
  double cosijk,cosTerm;
  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;

  static int firsttime=1;
  struct bondList *blist;
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  getBondList(&blist);
  

 
   
  if(firsttime){
    firsttime=0;
  }
  
  totPot=0.0;
  /*calculate forces*/
  for(ips=0;ips<atomsInPartSkin;ips++){
    i=partSkin[ips];
    bondListPosj=blist->head[i]; /*where the ngbrs start*/
    while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=blist->list[bondListPosj]; /*j is the ngbr*/
      bondListPosj++;
		
      drx_ij=pos.x[j]-pos.x[i]; 
      dry_ij=pos.y[j]-pos.y[i]; 
      drz_ij=pos.z[j]-pos.z[i];	
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		  
      if(rij2>minhBox2){
	PERIODIC(drx_ij,dry_ij,drz_ij);
	rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      }
      rij=sqrt(rij2);
      /*first the pair part*/
      /*here we sum over  i=0...N-1, j>i*/
      
      if(j>i){ /*dont doublecount!! in pair part*/
		  poti=KSTRETCH(i,j)*POW2(R0(i,j)-rij);
	totPot+=poti;
	pos.enPot[i]=0.5*poti;
	pos.enPot[j]=0.5*poti;
			
	if(calcForces){ /*j=n*/
	  forcej=2.0*KSTRETCH(i,j)*(R0(i,j)-rij)/rij;
	  /*xyza is acceleration, NOT force*/
	  pos.xa[i] -= forcej*drx_ij;
	  pos.ya[i] -= forcej*dry_ij;
	  pos.za[i] -= forcej*drz_ij;
			  
	  pos.xa[j] += forcej*drx_ij;
	  pos.ya[j] += forcej*dry_ij;
	  pos.za[j] += forcej*drz_ij;
	}
      }
		 
		 
      /*then the three body part */
      /*here we sum over i=0...N-1, j!=i, k!=i && k>j */
		 
      bondListPosk=bondListPosj; /*bondListPosj was increased already erlier so we dont need +1*/
      while(blist->list[bondListPosk]!=-1){ /*get second ngbr k>j*/
	k=blist->list[bondListPosk]; /*k is the ngbr*/
	bondListPosk++;
		  
	drx_ik=pos.x[k]-pos.x[i]; 
	dry_ik=pos.y[k]-pos.y[i]; 
	drz_ik=pos.z[k]-pos.z[i];	
	rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
		  
	if(rik2>minhBox2){
	  PERIODIC(drx_ik,dry_ik,drz_ik);
	  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
	}
	rik=sqrt(rik2);
		
		
	cosijk=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/(rik*rij);
	cosTerm=(cosijk-COS0(i));
	poti=KBEND(j,i,k)*POW2(cosTerm);
	pos.enPot[i]=poti;
	totPot+=poti;

		  
	if(calcForces){
	  varTemp1=-2.0*KBEND(j,i,k)*cosTerm;
	  varTemp2=varTemp1/(rik*rij);
	  varTemp3=varTemp1*cosijk/rij2;
	  varTemp4=varTemp1*cosijk/rik2;
			

	  forcej=varTemp2*drx_ik-varTemp3*drx_ij;
	  forcek=varTemp2*drx_ij-varTemp4*drx_ik;
	  pos.xa[j] += forcej;
	  pos.xa[k] += forcek;
	  pos.xa[i] += -forcej-forcek;
			 
	  forcej=varTemp2*dry_ik-varTemp3*dry_ij;
	  forcek=varTemp2*dry_ij-varTemp4*dry_ik;
	  pos.ya[j] += forcej;
	  pos.ya[k] += forcek;
	  pos.ya[i] += -forcej-forcek;
			 
	  forcej=varTemp2*drz_ik-varTemp3*drz_ij;
	  forcek=varTemp2*drz_ij-varTemp4*drz_ik;
	  pos.za[j] += forcej;
	  pos.za[k] += forcek;
	  pos.za[i] += -forcej-forcek;
	}
      }
    }
  }
  return totPot;
}
void putInteractingInTable(int atom,struct parameters *par,struct systemPos *pos,struct bondList *blist,int *table){
  int j,bondListPosj,ngbrListPosj;
  double drx_ij,dry_ij,drz_ij;
  double rij2;

  table[atom]=1;
  bondListPosj=blist->head[atom]; /*where the ngbrs start*/
  
  while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
    j=blist->list[bondListPosj]; /*j is the ngbr*/
    bondListPosj++;
    table[j]=1;
  }
  
  if(vff.includeRepPot){
    ngbrListPosj=nd.head[atom]; /*where the ngbrs start*/
    while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=nd.list[ngbrListPosj++]; /*j is the ngbr*/
      if(!table[j]){
	drx_ij=pos->x[j]-pos->x[atom]; 
	dry_ij=pos->y[j]-pos->y[atom]; 
	drz_ij=pos->z[j]-pos->z[atom];	
	PERIODIC(drx_ij,dry_ij,drz_ij);
	rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
	if( rij2 <vff.repR02 )
	  table[j]=1;
      }
    }
  }
}

void updateParts(struct parameters *par,struct systemPos pos,int *part,int *partSkin, int *atomsInPart,int *atomsInPartSkin,double updateTreshold){
  static int firsttime=1;
  int i,ip,j;
  struct bondList *blist;
  static int *table;
  int ngbrListPosj,bondListPosj;
  double drx,dry,drz;
  double r2;
  
  
  double minhBox2=MIN3(par->hBox[0],par->hBox[1],par->hBox[2]);
  minhBox2*=minhBox2; /*make sure its the ^2 .*/
  
  if(firsttime){
    table=calloc(pos.nAtoms,sizeof(int)); /* in this table we put a 1 if it should be in partSkin (the new one) */
    firsttime=0;
  }
  getBondList(&blist);
  
  /* first new part has to be created */
  
  memset(table,0,pos.nAtoms*sizeof(int));
  
  for(ip=0;ip<(*atomsInPart);ip++){
    i=part[ip];
    drx=pos.x[i]-pos.prevx[i]; 
    dry=pos.y[i]-pos.prevy[i]; 
    drz=pos.z[i]-pos.prevz[i];	
    r2=LENGTH2(drx,dry,drz);
    if(r2>minhBox2){
      PERIODIC(drx,dry,drz);
      r2=LENGTH2(drx,dry,drz);
    }
    
    if( r2 >updateTreshold )
      putInteractingInTable(i,par,&pos,blist,table);
    else
      table[i]=1; /*it has not moved, lets not include its neighbours*/
  }
  
  /* now write this to the part list */
  j=0;
  for(i=0;i<pos.nAtoms;i++)
    if(table[i])
      part[j++]=i;
  (*atomsInPart)=j;
  
  

  /* then new part+skin has to be created */
  memset(table,0,pos.nAtoms*sizeof(int));
  /* first put all bonded interacting particle */
  for(ip=0;ip<(*atomsInPart);ip++)
    putInteractingInTable(part[ip],par,&pos,blist,table);
    
  /* now write this to the partSkin table */
  j=0;
  for(i=0;i<pos.nAtoms;i++)
    if(table[i])
      partSkin[j++]=i;
  
  (*atomsInPartSkin)=j;
  
}
void reinitAfterStep(struct parameters *par,struct systemPos pos,int stepAccepted){
  if(stepAccepted){
    copyNgbrs(&ndLastStep,&nd,pos.nAtoms);
  }
  else
    copyNgbrs(&nd,&ndLastStep,pos.nAtoms);
}

double partForcesVff(struct parameters *par,struct systemPos pos,int *part, int *atomsInPart,int calcForces,int updatePart,double updateTreshold){
  static int firsttime=1;
  static int *partSkin;
  static int atomsInPartSkin;

  double totPot=0.0;
  
  if(firsttime){
    partSkin=malloc(sizeof(int)*pos.nAtoms);
    atomsInPartSkin=0;
    firsttime=0;
  }
  
  
  if(calcForces){ /*this is not good if we want to have more than one pot at the same time!!!*/
    memset(pos.xa,0,pos.nAtoms*sizeof(double));
    memset(pos.ya,0,pos.nAtoms*sizeof(double));
    memset(pos.za,0,pos.nAtoms*sizeof(double));
  }

  if(updatePart==2){				  /* new minimization started!!, the part arrays have to be initialized */
    memcpy(partSkin,part,(*atomsInPart)*sizeof(int)); /* when the parts are updated part will become partskin and partskin will be updates..  */
    atomsInPartSkin=(*atomsInPart);
  
  }

  if(updatePart)
    updateParts(par,pos,part,partSkin,atomsInPart,(&atomsInPartSkin),updateTreshold);
  
  /*calculate forces*/
  
  if(vff.includeVffPot)  /*dont know why one would want to exclude this part but...*/
    //  totPot+=partKeatingNotable(par,pos,partSkin,atomsInPartSkin,calcForces); 
  	 totPot+=partKeating(par,pos,partSkin,atomsInPartSkin,calcForces);
  
  if(vff.includeRepPot)
    totPot+=partRepulsive(par,pos,partSkin,atomsInPartSkin,calcForces);
  
  if(vff.includeSuboxidePot)
    totPot+=partSuboxide(par,pos,part,*atomsInPart);
  
  if(vff.includeDanglbondPot)
    totPot+=partDanglBonds(par,pos,part,*atomsInPart);
  
  return totPot;
  

}














