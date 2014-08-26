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



#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "shared.h"
#include "mcSubs.h"
#include "initialize.h"
#include "fileio.h"
#include "parallel.h"
#include "bondList.h"
#include "unistd.h"
#include "time.h"


/* This module contains functions that make file input and output easier
   the structure fileIo contains information about the names of the files

   struct fileIo
   char inFileName[256];    the name of the file from where we read in the system
   char inFileFormat[256];  its type
   char outFileName[256];   the name of the file where we write in the system
   char outFileFormat[256]; its type
   char comment[1024];      a comment of the system 
   char opendxFileName[256];   the name of the file where we write fdata for the opendxclient program
   char opendxFileFormat[256]; its type
   char outputFileName[256];   the name of the file where we write some information at each step, like temp energy
	
   the functions which are defined in the interface are
	
   void initStructFileIo(struct fileIo *fio);
   void openSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par); 
   void writeSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par);
   void writeOpendxData(struct fileIo *fio,struct systemPos *pos,struct parameters *par);
   void writeOutput(struct fileIo *fio,char *outString);
   void openParameterFile(char *filename);
   void closeParameterFile();
   int getParValue(char *parameter,void *value,char *formatString);
	
   initStructFileIo:
   puts some default values int the fields of the structure

   openSystem:
   reads in the system in file inFileName of format inFileFormat into structs pos and par
	
   writeSystem:
   writes the system to the file outFileName of format inFileFormat from the structs pos and par
	
   writeOpendxData:
   writes the system the file opendxFileName of format  opendxFileFormat from the structs pos and par
	
   writeOutput:
   writes th given string outString the file outputFileName. The first time this function is called the
   file is opened automatically
	
   The last three functions are used to handle parameterfiles. Parameterfiles are files
   where we have a lines with parameterlabels and their value. Comments can be put in
   with // . The functions are quite forgiving and it oesnt mather if there are empty lines, lines with normal
   text and os on.
   example:

   //alphas value is 0.44
   alpha 0.44

   To read parameterfiles we first have to open the file with openParameterFile(char *filename)
   after this we can get values with getParValue(char *parameter,void *value,char *formatString)
   to read aplha we would give the command: 
	
   int i;
   double a;
   i=getParValue("alpha",&a,"%lf");
	
   the values of i:  0  The value was read correctly
   1  The parameter was not found in the file
   2  Some other error (ex the value was not given but the label was)

   Only one parameter file can be open at a time so it has to be closed with  void closeParameterFile();
   after all the parameters have been read
	
*/

/*internal structs and defines*/
#define PARAMSIZE 256
struct paramVal{
  char parameter[PARAMSIZE];
  char value[PARAMSIZE];
};

/*internal functions*/

void writePdb(struct systemPos pos,struct parameters *par,char fname[],char comment[]);
void openXyz(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]);
void openPdb(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]);
void writeXyz(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]);


/*internal globals*/
static struct paramVal *pv=NULL;
static char paramFileName[256];

/*==============FUNCTIONS==================*/

void initStructFileIo(struct fileIo *fio){
  strcpy(fio->inFileName,"in.xyz");
  strcpy(fio->inFileFormat,"xyz");
  strcpy(fio->outFileName,"out.xyz");
  strcpy(fio->outFileFormat,"xyz");
  strcpy(fio->comment,""); 
}


void writeSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par, char *comment){
  int rank;
  static int callNum=0;			  /* the number of times this function has been called, used to give new names to the outfile */
  char fileName[256];
  char buffer[200];
  int isGzipped;
  rank=getRank(); 
  sprintf(fileName,"%03d_%s",callNum,fio->outFileName);

  if(strstr(fileName,".gz")==NULL)
    isGzipped=0;
  else{
    isGzipped=1;
    *(strstr(fileName,".gz"))='\0'; /*  take away gz from filename, the filewriting function s first have to write unzipped versions of the files*/
  }
  
  if(rank==0){
    if(strncmp(fio->outFileFormat,"xyz",3)==0)
      writeXyz(pos,par,fileName,fio->outFileFormat,comment);
    else if(strncmp(fio->outFileFormat,"pdb",3)==0)
      writePdb(*pos,par,fileName,comment);
	
    else
      printf("unknown fileformat, could not write system\n");
  }
  	 	
  if(isGzipped){
    sprintf(buffer,"gzip -f %s",fileName);
    system(buffer);			  /* now zip the file.. */
	 
  }
  
  callNum++;

}




void openSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par){
  int rank;
  int isGzipped;
  char filename[200];
  char buffer[200];
  rank=getRank();
  strcpy(filename,fio->inFileName);
  
  if(rank==0) printf("Reading system from file %s...",filename);fflush(stdout);
  
  if(rank==0){  
    if(strstr(filename,".gz")==NULL)
      isGzipped=0;
    else{
      isGzipped=1;
      sprintf(buffer,"gunzip %s",filename);
      system(buffer);			  /* unzrip packed file */
      *(strstr(filename,".gz"))='\0'; /*  take away gz from filename */
    }
	 

    if(strncmp(fio->inFileFormat,"xyz",3)==0)
      openXyz(pos,par,filename,fio->inFileFormat,fio->comment); 
    else if(strncmp(fio->inFileFormat,"pdb",3)==0)
      openPdb(pos,par,filename,fio->inFileFormat,fio->comment);
    else
      printf("unknown fileformat (%s), could not read  system\n",fio->inFileFormat);
	 	
    if(isGzipped){
      sprintf(buffer,"gzip %s",filename);
      system(buffer);			  /* zip file again file */
    }	 
  }

#ifdef PARALLEL
  MPI_Bcast(&pos->nAtoms,1,MPI_INT,0,MPI_COMM_WORLD);
  if(rank!=0) initStructPos(pos,pos->nAtoms);  
  MPI_Bcast(pos->x,3*pos->nAtoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(pos->aType,pos->nAtoms,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(pos->aFlag,pos->nAtoms,MPI_LONG,0,MPI_COMM_WORLD); 
  MPI_Bcast(&(par->box[0]),3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  setSystemSize(par);
#endif
  
  if(rank==0)   printf(" ready. System has %d atoms, box=%g %g %g\n",pos->nAtoms,par->box[0],par->box[1],par->box[2]);
}


void openXyz(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]){
  
  char xyzatom[3];
  char buf[1024];
  int i,n,nAtoms;
  FILE *fp;
  
  /*open original file*/
  fp=fopen(fname,"r");
  if (fp==NULL){ 
    printf("Could not open %s\n",fname);
    exit(0);
  }
  fgets(buf,256,fp);
  sscanf(buf,"%d",&nAtoms);
  initStructPos(pos,nAtoms);  
  fgets(comment,1024,fp);
  for (i=0;i<pos->nAtoms;i++) {
    pos->aFlag[i]=0;
  }

  if(strcmp(ftype,"xyz")==0) /*normal xyz*/
    for (i=0;i<pos->nAtoms;i++) {
      fgets(buf,256,fp);
      n=sscanf(buf,"%s %lg %lg %lg %d",xyzatom,pos->x+i,pos->y+i,pos->z+i,pos->aType+i);
    }
  

  else if (strcmp(ftype,"xyz_flag")==0)
    for (i=0;i<pos->nAtoms;i++) {
      fgets(buf,256,fp);
      n=sscanf(buf,"%s %lg %lg %lg %d %ld",xyzatom,pos->x+i,pos->y+i,pos->z+i,pos->aType+i,pos->aFlag+i);
    } 
  else if (strcmp(ftype,"xyz_noatype_sio")==0)
    for (i=0;i<pos->nAtoms;i++) {
      fgets(buf,256,fp);
      n=sscanf(buf,"%s %lg %lg %lg ",xyzatom,pos->x+i,pos->y+i,pos->z+i);
      if(strcmp("Si",xyzatom))
		  pos->aType[i]=1;
      else
		  pos->aType[i]=0;
		
    }
  else{
    printf("unkown file format\n");
    exit(0);
  }
  fclose(fp);
}


void openPdb(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]){
  
  double bx,by,bz;
  char buf[1024];
  char line[10];
  char dummy[20];
  char xyzAtom[3];
  int i,n,nAtoms;
  FILE *fp;
  
  for (i=0;i<pos->nAtoms;i++) 
    pos->aFlag[i]=0;
  
  /*open original file*/
  fp=fopen(fname,"r");
  if (fp==NULL){ 
    printf("Could not open %s\n",fname);
    exit(0);
  }

  fgets(buf,256,fp);
  sscanf(buf,"%s %s %d",line,dummy,&nAtoms);
  if(strcmp(dummy,"NATOMS")){
    if(getRank()==0) printf("could not read num of atoms from pdb file\n");
    exit(0);
  }
  initStructPos(pos,nAtoms);  
  
  fgets(buf,256,fp);
  sscanf(buf,"%s %s %lg %lg %lg",line,dummy,&bx,&by,&bz);
  if(!strcmp(dummy,"BOX")){	  /* if the box size was defined  */
    if(getRank()==0)
      printf("Box size read from pdbfile, is %g %g %g\n",bx,by,bz);

    par->box[0]=bx;
    par->box[1]=by;
    par->box[2]=bz;
    setSystemSize(par);
  }
  

  
  
  do{			 /* go to start of atoms */
    fgets(buf,256,fp);
    sscanf(buf,"%s",line);
  }  while(strcmp(line,"ATOM"));
  

  
  for (i=0;i<pos->nAtoms;i++) {
    n=sscanf(buf,"%*s %*d %s %*s %*d %lg %lg %lg %*g %*g %ld",xyzAtom,pos->x+i,pos->y+i,pos->z+i,pos->aFlag+i);
   
    if(n==4)
      pos->aFlag[i]=0; /*aFlag not part of pdb standard, might not be there..*/
    
    if(strcmp(xyzAtom,"Si")==0) /* abit of kludge */
      pos->aType[i]=SI_ATYPE;
    else if(strcmp(xyzAtom,"O")==0)
      pos->aType[i]=O_ATYPE;
    fgets(buf,256,fp);
  }
  
  fclose(fp);
}



void writePdb(struct systemPos pos,struct parameters *par,char fname[], char comment[]){
  FILE *fp;
  FILE *fp2;
  int i,j,ngbrListPosj;
  int *head;
  int *list;
  char host[64];
  char bondfname[200];
  double drx_ij,dry_ij,drz_ij;
  double r2;
  int numOfBonds;
  time_t now;
  struct bondList *blist;
  getBondList(&blist);
  getBondInfo(&head,&list);	  /* ugly, we could just use blist which contains these, oh well */
	
  gethostname(&(host[0]),(int)sizeof(host));
  now=time(NULL);
  strncpy(bondfname,fname,190);
  strcat(bondfname,".bonds");
  printf("Writing out %d atom coordinates to %s, bonds to %s...",pos.nAtoms,fname,bondfname);fflush(stdout);
  fp=fopen(fname,"wt");
  
  
  fprintf(fp,"REMARK     NATOMS %d\n",pos.nAtoms);
  fprintf(fp,"REMARK     BOX    %g %g %g \n",par->box[0],par->box[1],par->box[2]);
  fprintf(fp,"REMARK     HOST   %s\n",host);
  fprintf(fp,"REMARK     PWD    %s\n",getenv("PWD"));
  fprintf(fp,"REMARK     USER   %s\n",getenv("USER"));
  fprintf(fp,"REMARK     DATE   %s",ctime(&now));
  fprintf(fp,"REMARK     PROG   bomc (%s %s)\n",__DATE__,__TIME__);
  fprintf(fp,"REMARK     %s\n",comment);
  
  /*SYstem info*/
  fprintf(fp,"USERSY %d %g %g %g\n",pos.nAtoms,par->box[0],par->box[1],par->box[2]);
  /*SImulation info */
  /*
    fprintf(fp,"USERSI %d %g %g %g",pos.nAtoms,par->box[0],par->box[1],par->box[2]);
  */

  for(i=0;i<pos.nAtoms;i++){ 
	  	  
    numOfBonds=0;
    ngbrListPosj=head[i]; /*where the ngbrs start*/
    while(list[ngbrListPosj++]!=-1) /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      numOfBonds++;
	  
    fprintf(fp,"ATOM  ");
    fprintf(fp,"%5d ",i);
    fprintf(fp,"%4s",par->aName[pos.aType[i]]);
    fprintf(fp," MOL     1    ");
    fprintf(fp,"%8.3f",pos.x[i]);
    fprintf(fp,"%8.3f",pos.y[i]);
    fprintf(fp,"%8.3f",pos.z[i]);
    fprintf(fp,"%6.2f",(double)numOfBonds);
    fprintf(fp,"%6.2f",pos.enPot[i]*TOEV);
    fprintf(fp,"%5ld\n",pos.aFlag[i]); /*non standard*/
	
	  	  
  }
  
  fprintf(fp,  "TER   %5d      MOL     1\n",pos.nAtoms);
  
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    ngbrListPosj=head[i]; /*where the ngbrs start*/
    while(list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=list[ngbrListPosj]; /*j is the ngbr*/
      ngbrListPosj++;	
      if(i<j){ /*dont write unneccessary bons...*/
	drx_ij=pos.x[j]-pos.x[i]; 
	dry_ij=pos.y[j]-pos.y[i]; 
	drz_ij=pos.z[j]-pos.z[i];	
	r2=LENGTH2(drx_ij,dry_ij,drz_ij);
	PERIODIC(drx_ij,dry_ij,drz_ij);
	if(fabs(r2-LENGTH2(drx_ij,dry_ij,drz_ij))<0.01 ) /*we want to exclude atoms that are in different parts of the system*/
	  fprintf(fp,  "CONECT%5d     %5d\n",i,j);
      }
    }
  }
  

 for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i, write additional AtomicInfo*/
    fprintf(fp,"USERAI%5d  %g  \n",i,pos.enPot[i]*TOEV);
  }

  fp2=fopen(bondfname,"wt");
  for(i=0;i< pos.nAtoms;i++){ /*loop over atoms i*/
    /*loop ov er atoms i, write bonding information (all of it!) BondInfo*/
    numOfBonds=0;
    ngbrListPosj=head[i]; /*where the ngbrs start*/
    while(list[ngbrListPosj++]!=-1) /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      numOfBonds++;
    fprintf(fp2,"%d ",numOfBonds);
    fprintf(fp,"USERBI%5d  %d %d ",i,numOfBonds,blist->danglBonds[i]);
    ngbrListPosj=head[i]; /*where the ngbrs start*/
    while(list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      fprintf(fp2,"%d ",list[ngbrListPosj]); 
      fprintf(fp,"%d ",list[ngbrListPosj++]); 
    }
    fprintf(fp2,"%d \n",blist->danglBonds[i]);
    fprintf(fp,"\n");
  }

  fclose(fp2);

  fprintf(fp,"END\n");
  fclose(fp);
  printf("ready\n");
}


void writeXyz(struct systemPos *pos,struct parameters *par,char fname[],char ftype[],char comment[]){
  FILE *fp;
  int i;
	
  printf("Writing out %d atom coordinates to %s...",pos->nAtoms,fname);fflush(stdout);
  fp=fopen(fname,"wt");
  fprintf(fp,"%d\n",pos->nAtoms);
  fprintf(fp,"%s boxsize %g %g %g periodic \n",comment,par->box[0],par->box[1],par->box[2],par->periodic[0],par->periodic[1],par->periodic[2]);

  if(strcmp(ftype,"xyz")==0) /*normal xyz*/
    for(i=0;i<pos->nAtoms;i++) 
      fprintf(fp,"%s %16.11f %16.11f %16.11f %d\n",par->aName[pos->aType[i]],pos->x[i],pos->y[i],pos->z[i],pos->aType[i]);
  else if (strcmp(ftype,"xyz_flag")==0)
    for(i=0;i<pos->nAtoms;i++) 
      fprintf(fp,"%s %16.11f %16.11f %16.11f %d %ld\n",par->aName[pos->aType[i]],pos->x[i],pos->y[i],pos->z[i],pos->aType[i],pos->aFlag[i]);
  else
    printf("unkown file format,nothing written\n");

  fclose(fp);
  printf("ready\n");
}


/* the  following three functions amkes it  possible to read parameters from a file*/

void openParameterFile(char *filename){
  int i;
  char line[256];
  int numLines; /*this is the number of lines with values */
  FILE *fp;
  int rank;
  rank=getRank();
  if(rank==0){
    if(pv!=NULL){
      printf("another paramter file(%s) is already open, one at the time....\n",paramFileName);
      exit(0);
    }
    /*lets first calculate how many lines the file contains*/
    /*this is a bit of a kludge, first we open file and get the number of lines, the we open it again and read the
      items..*/
	 
    strcpy(paramFileName,filename);
    fp=fopen(filename,"r");
    if (fp==NULL){ 
      printf("Could not open %s\n",filename);
      exit(0);
    }
    numLines=0;
    while(fgets(line,256,fp)!=NULL)
      if(strncmp(line,"//",2)!=0 && strcmp(line,"")!=0) /*checks if we have empty or comment lines*/ 
	numLines++;
    fclose(fp);
	 
    pv=malloc(sizeof(struct paramVal)*(numLines+1));
    /*lets read in the values*/
    fp=fopen(filename,"r");  
    i=0;
    while(fgets(line,256,fp)!=NULL)
      if(strncmp(line,"//",2)!=0 && strcmp(line,"")!=0){ /*checks if we have empty or comment lines*/ 
	sscanf(line,"%s %s",pv[i].parameter,pv[i].value);
	i++;
      }
    fclose(fp);
    strcpy(pv[i].parameter,"theend"); /*is used to detect the end of the array pv */
  }
  
#ifdef PARALLEL
  MPI_Bcast(&numLines,1,MPI_INT,0,MPI_COMM_WORLD);
  if(rank!=0)
    pv=malloc(sizeof(struct paramVal)*(numLines+1));
  for(i=0;i<numLines+1;i++){
    MPI_Bcast(pv[i].parameter,PARAMSIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(pv[i].value    ,PARAMSIZE,MPI_CHAR,0,MPI_COMM_WORLD);
  } 
#endif

}

void closeParameterFile(void){
  free(pv);
  pv=NULL;
}

int getParValue(char *parameter,void *value,char *formatString){
  int i=0,j;
  
  while(strcmp(pv[i].parameter,"theend")!=0){
    if(strcmp(pv[i].parameter,parameter)==0){
      j=sscanf(pv[i].value,formatString,value);
      if(j==1)
	return 0; /*we found the value*/
      else{
	printf("could not read value of parameter %s. Should be of type %s\n",parameter,formatString);
	return 2; /*could not get the value*/
      }
    }
    else
      i++;
  }
  printf("could not find parameter %s\n",parameter);
  return 1; /*we did not find the parameter in the list*/
}


void readBondFile(struct fileIo *fio,struct systemPos *pos,struct parameters *par){
  struct bondList *blist;
  static FILE *fp;
  char buf[256];
  int i,j,n;
  int bond[4];
  int numOfBonds;
  int totNumOfBonds;
  int danglBonds;
  int bondListPosj;
  
  getBondList(&blist);
  
  /*copy bonds from ndBlist, except O-O bonds*/
  if(getRank()==0){				  /* if rank is 0 or this is a serial prog */
    printf("reading bondlist...");
    fflush(stdout);
    
    fp=fopen(fio->inBondFile,"r"); 
    totNumOfBonds=0;
    for(i=0;i< pos->nAtoms;i++){ /*loop over atoms i*/
      bondListPosj=blist->head[i];
		
      fgets(buf,256,fp);
      n=sscanf(buf,"%d %d %d %d %d %d ",&numOfBonds,&bond[0],&bond[1],&bond[2],&bond[3],&danglBonds);
      totNumOfBonds+=numOfBonds;
		
      if(numOfBonds!=4) /* if less than 4 bonds, dangl bond info went ot the wrong place */
	danglBonds=bond[numOfBonds]; 
		
      blist->danglBonds[i]=danglBonds;
		
      for(j=numOfBonds;j<4;j++)
	bond[j]=-1;	/* no bonds */

      for(j=0;j<4;j++) /* lets finally write info tp bondlist */
	blist->list[bondListPosj++]=bond[j];
      blist->list[bondListPosj++]=-1; /* end marker */
    }

    fclose(fp);
    printf("ready, read in  %d bonds %g per atom\n",totNumOfBonds,(double)totNumOfBonds/pos->nAtoms);
  }
  
#ifdef PARALLEL
  MPI_Bcast(blist->head,blist->arraySize,MPI_INT,0,MPI_COMM_WORLD);
#endif
}
