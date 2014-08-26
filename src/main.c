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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include "shared.h"
#include "fileio.h"
#include "potentialWww.h"
#include "initialize.h"
#include "createLattice.h"
#include "mcSubs.h"
#include "miscSubs.h"
#include "conjgrad.h"
#include "parallel.h"
#include "randFunc.h"
#include "bondList.h"                     
#include "rings.h"        
#include "wwwMc.h"


/*

Sebastian von Alfthan  2003

Bomc: program for calculating WWW MC simulations


*/

/*internal functions*/
void readParameterFile(int argc,char *argv[]);
void programAborted(int signal);




/*globals in module*/
static struct systemPos pos;           /*contains information about each atom such as 
					 position,speed,force,charge,type and so on,  defined in shared.h*/
static struct parameters par;          /*contains parameters which describe the system, such as boxsize
					 periodicity,name of a certain atomtype, mass of a certain atomtype
					 and so on, defined in shared.h*/
static struct fileIo fio;              /*contains information about the input,output of the program sunch as name of
				          files, fileetypes and so on.  defined in fileio.h*/
static struct simuPar simuPar;            /*contains information about how to perform the simulation. defined in main.c*/




int main(int argc,char *argv[]){
  time_t startTime;				  /*  when the simulation (anneal part) started */
  char buffer[512];
  int rank;


#ifdef PARALLEL
  /*init mpi calculation*/
  MPI_Init(&argc,&argv); 
#endif
  /* catch signals, if caught run programAborted that wirtes out the current situation of system, except TERM which is ignored */
  
  signal(SIGABRT,programAborted);
  signal(SIGFPE,programAborted);
  signal(SIGILL,programAborted);
  signal(SIGINT,programAborted);
  signal(SIGSEGV,programAborted);
  signal(SIGHUP,SIG_IGN);		  /* ignore hangup */
 
  rank=getRank();
  if(getRank()==0)  printf("bomc compiled %s %s\n", __DATE__,__TIME__);

  startTime=time(NULL);
  
  
  
  /*initialize the structures (puts some default values in them)*/
  initStructParameters(&par);
  initStructFileIo(&fio);
  
  /*read the parameter file (and initialize a bunch of stuff)*/
  readParameterFile(argc,argv);
  
  /*init bonds*/
  initBondList(pos);
  if(simuPar.createBonds)
    connectAtoms(pos,&par,simuPar.bondCreationR0);
  else{
    readBondFile(&fio,&pos,&par);
    deleteLongBonds(pos,&par,5.0); /*
				    * delete bonds longer than 5.0 A, unphysical 
				    * (useful if one wants to transform a periodic 
				    * system to nonperiodic)
				    */
    
  }
  initwwwPot(&par,pos);

  
  
    /*
	 check force to see if correct, that is the force is the derivative of the potential
  */

  if(simuPar.checkBulkmodulus){
	 calculateBulkModulus(pos,&par);
	 exit(0);
  }
  if(simuPar.checkPotential){
	 checkForceCalc(pos,&par);
    exit(0);
  }
  
  sprintf(buffer,"original input");
  writeSystem(&fio,&pos,&par,buffer);	  
  prepareInitialSystem(pos,&par,&simuPar);
  sprintf(buffer,"Initial system used in simulations");
  writeSystem(&fio,&pos,&par,buffer);	   
  
 
  /*this should make sure the different processors start from diff seeds, in the mc part we want them to behave differently*/  
  /*simulation ready*/
  initRandNum(3,*par.seed+rank);
  mcSimulation(pos,&par,&simuPar,&fio);

  if(getRank()==0){
    int days,hours,minutes,seconds;
    double potVal;
    potVal=wwwPot(&par,pos,0);
    convertSeconds((int)difftime(time(NULL),startTime),&days,&hours,&minutes,&seconds);
    sprintf(buffer,"Simulation ready  kT:g%eV pot:%geV executiontime:%dd %dh %dm %ds"
	    ,par.kT*TOEV,potVal*TOEV/pos.nAtoms,days,hours,minutes,seconds);  
    if(getRank()==0) printf("%s\n",buffer);
    writeSystem(&fio,&pos,&par,buffer);
    if(simuPar.calculateDynMatrix)
      calculateDynamicalMatrix(pos,&par);
  }
  
  
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}


/*reads the parameter file*/
void readParameterFile(int argc,char *argv[]){
  int haveSys=0;
  int error=0;
  char input[20];
  double temp;
  int  rank=getRank();
  /*lets get the name of the parameter file and open it*/
  
  openParameterFile("main.par");
  error+=getParValue("seed",par.seed,"%li");
  initRandNum(3,*par.seed-1);/*this should make sure the different processors start from same seed, teh randomization should be the same in all processors*/
  
  /*read input method of system, can be file or create*/
  error+=getParValue("input",input,"%s");
  
  /*we get the sys from file*/
  if(strcmp(input,"file")==0){ 
    error+=getParValue("infile_name",fio.inFileName,"%s");
    error+=getParValue("infile_format",fio.inFileFormat,"%s");
    error+=getParValue("box_x",&par.box[0],"%lf");
    error+=getParValue("box_y",&par.box[1],"%lf");
    error+=getParValue("box_z",&par.box[2],"%lf");  
    setSystemSize(&par);
    
    if(error==0){
      openSystem(&fio,&pos,&par);
      haveSys=1;
    }
    else{
      if(rank==0)	printf("something wrong with the values for infile!\n");
      exit(0);
    }
  }
  
  /*we have to create the sys*/
  else if(strcmp(input,"create")==0){
    double density;
    double boxx,boxy,boxz;
    
    char sys[20];
	 /*read in which system we want to create and the size of it in unitcells*/
    error+=getParValue("sys_name",sys,"%s");
    error+=getParValue("density",&density,"%lf");
    error+=getParValue("create_box_x",&boxx,"%lf");
    error+=getParValue("create_box_y",&boxy,"%lf");
    error+=getParValue("create_box_z",&boxz,"%lf");
    
    if(strcmp(sys,"Si")==0){
      createZnBLattice(boxx,boxy,boxz,1,1,density,&pos,&par);
      haveSys=1;
	 } 
    if(strcmp(sys,"SiO2")==0){
      createIbcLattice(boxx,boxy,boxz,1,0,2,density,&pos,&par); /* create ideal beta cristobalite */
      haveSys=1;
    } 
    if(strcmp(sys,"SiO")==0){
      createIbcLattice(boxx,boxy,boxz,1,0,1,density,&pos,&par); 
      haveSys=1;
    } 

    if(strcmp(sys,"Sirand")==0){
      createRandLattice(boxx,boxy,boxz,density,1,1,1,&pos,&par);
      haveSys=1;
    }
    if(strcmp(sys,"SiO2rand")==0){
      createRandLattice(boxx,boxy,boxz,density,0,1,2,&pos,&par);
      haveSys=1;
	 }
    if(strcmp(sys,"SiOrand")==0){
      createRandLattice(boxx,boxy,boxz,density,0,1,1,&pos,&par); /* what is correct density for SiO ?? */
      haveSys=1;
    }
  }

  
  error+=getParValue("createBonds",&simuPar.createBonds,"%d");  
  error+=getParValue("inBondFile",fio.inBondFile,"%s");
  error+=getParValue("bondCreationR0",&simuPar.bondCreationR0,"%lf");  


  if(haveSys==0 || error!=0){
	 	if(rank==0)printf("something wrong with the system input paramters!!\n");
	 exit(0);
  } 
  
  error+=getParValue("outfile_name",fio.outFileName,"%s"); 
  error+=getParValue("outfile_format",fio.outFileFormat,"%s");
  
  if(error){
	 if(getRank()==0) printf("give outfile name and format\n");
	 exit(0);
  }
  
    
  error+=getParValue("result_period",&temp,"%lf"); 
  simuPar.resultPeriod=temp*pos.nAtoms;

  error+=getParValue("periodic_x",&par.periodic[0],"%d");
  error+=getParValue("periodic_y",&par.periodic[1],"%d");
  error+=getParValue("periodic_z",&par.periodic[2],"%d"); 

  error+=getParValue("varVol_x",&par.varVol[0],"%d");
  error+=getParValue("varVol_y",&par.varVol[1],"%d");
  error+=getParValue("varVol_z",&par.varVol[2],"%d"); 
		 	 
  
  error+=getParValue("rand_iter",&temp,"%lf "); 
  simuPar.randIter=(int)(temp*pos.nAtoms);

 			
  error+=getParValue("pressure",&par.pressure,"%lf");
  par.pressure*=PASCAL;
 
  error+=getParValue("oxygenDiffuseProb",&simuPar.oxDiffuseProb,"%lf"); 
  error+=getParValue("bondDiffuseProb",&simuPar.bondDiffuseProb,"%lf"); 
  error+=getParValue("bondSwitchProb",&simuPar.bondSwitchProb,"%lf"); 
  error+=getParValue("bondCreateProb",&simuPar.bondCreateProb,"%lf"); 
  error+=getParValue("bondBreakProb",&simuPar.bondBreakProb,"%lf"); 
  
  /* lets normalize these */
  temp=simuPar.oxDiffuseProb+simuPar.bondDiffuseProb+simuPar.bondSwitchProb+simuPar.bondCreateProb+simuPar.bondBreakProb;
  simuPar.oxDiffuseProb/=temp;
  simuPar.bondBreakProb/=temp;
  simuPar.bondSwitchProb/=temp;
  simuPar.bondCreateProb/=temp;
  simuPar.bondDiffuseProb/=temp;

  
  error+=getParValue("volstep_prob",&simuPar.volProb,"%lf"); 
  error+=getParValue("volstep_start_kT",&simuPar.volStartkT,"%lf"); 
  simuPar.volStartkT*=EV;

  /*  error+=getParValue("useCg",&simuPar.useCg,"%d");*/
  
  
  error+=getParValue("skin",&par.skin,"%lf");

  error+=getParValue("checkBulkmodulus",&simuPar.checkBulkmodulus,"%d");
  error+=getParValue("checkPotential",&simuPar.checkPotential,"%d");
  error+=getParValue("calculateDynMatrix",&simuPar.calculateDynMatrix,"%d");
  
  if(error){
    if(rank==0)	 printf("there were errors when getting parameters!\n");
    exit(0);
  } 
  closeParameterFile();
  
}
void programAborted(int sig){
  char buffer[200];
  signal(sig,SIG_DFL);		  /* if going into writeSystem evokes again signal we do not want to go into an infinite loop  */
  sprintf(buffer,"Simulation aborted with signal %d.",sig);  
  if(getRank()==0) printf("%s\n",buffer);
  fflush(stdout);
  writeSystem(&fio,&pos,&par,buffer);
  /* now it wiMll return back to the program and probably crash and burn :) */
}
