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
#include "tabulatedValues.h"
#include "shared.h"
#include <stdio.h>
#include <math.h>

void readTable(char fname[],struct valueTable *valTable){
  FILE *fp;
  char buffer[512];
  int i;
  double x2,x1;

  fp=fopen(fname,"r");
  valTable->elements=0;
  while(fgets(buffer,512,fp)!=NULL)
	 valTable->elements++;
  fclose(fp);
  
  
  valTable->y=malloc(sizeof(double)*valTable->elements);
  
  fp=fopen(fname,"r");
  fgets(buffer,512,fp);
  sscanf(buffer,"%lg %lg",&x1,valTable->y);
  fgets(buffer,512,fp);
  sscanf(buffer,"%lg %lg",&x2,valTable->y+1);
  valTable->xdelta=x2-x1;
  valTable->ixdelta=1.0/valTable->xdelta;
  
  for(i=2;i<valTable->elements;i++){
	 fgets(buffer,512,fp);
	 sscanf(buffer,"%*lg %lg",valTable->y+i);
  }

  for(i=0;i<valTable->elements;i++)
	 valTable->y[i]*=EV;

  fclose(fp);
}



/*
Now done using linear approximation...
*/


double tabulatedVal(double x,struct valueTable valTable){
  int n,m;
  double val;
  n=x/valTable.xdelta;
  m=n+1;
  if(m<valTable.elements){
	 val= valTable.y[n]+(x-n*valTable.xdelta)/valTable.xdelta*(valTable.y[m]-valTable.y[n]);
	 /*	 printf("x:%g val:%g  n:%d m:%d valn:%g valm:%g \n",x,val*TOEV,n,m,valTable.y[n]*TOEV,valTable.y[m]*TOEV);*/
	 return valTable.y[n]+(x-n*valTable.xdelta)*valTable.ixdelta*(valTable.y[m]-valTable.y[n]);
  
  }
  else
	 return 0.0;
}

double tabulatedDeriv(double x,struct valueTable valTable){
  int n,m;
  n=x/valTable.xdelta;
  m=n+1;
  if(m<valTable.elements){
	 return (valTable.y[m]-valTable.y[n])*valTable.ixdelta;
  }
  else
	 return 0.0;

}

