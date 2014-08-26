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
#include <stdio.h>
#include <memory.h>
#include "memSubs.h"


double **allocMatrix2(int n,int m){
  double **a;
  int i;
  a=malloc(sizeof(double*)*n);
  (*a)=malloc(sizeof(double)*m*n);
  
  if( a[0]==NULL ||a==NULL){
	 printf("could not allocate memory\n");
	 exit(0);
  }
  for(i=1;i<n;i++)
	 a[i]=a[i-1]+m;
  return a;
}

double ***allocMatrix3(int n,int m,int l){
  double ***a;
  int i;
  
  a=malloc(sizeof(double**)*n);
  a[0]=malloc(sizeof(double *)*n*m);
  a[0][0]=malloc(sizeof(double)*n*m*l);
  if(a[0][0]==NULL || a[0]==NULL ||a==NULL){
	 printf("could not allocate memory\n");
	 exit(0);
  }
		
  for(i=1;i<n;i++)
	 a[i]=a[i-1]+m;
  for(i=1;i<n*m;i++)
	 a[0][i]=a[0][i-1]+l;
  return a;
}
