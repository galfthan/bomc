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
#include <stdlib.h>
#include <string.h>

typedef struct {
  long type;
  long mrows;
  long ncols;
  long imagf;
  long namelen;
} Fmatrix;


int writeMat(char *name,int rows,int cols,int isImag,double *realData,double *imagData){
  Fmatrix x;
  int mn;
  FILE *fp;
  char matrixName[64];                                                                                
  char fileName[68];                                                                                
  
  strncpy(matrixName,name,63);
  matrixName[63]='\0';			  /* in case length of name is bigger than 63.. */
  strcpy(fileName,matrixName);
  
  fp=fopen(strcat(fileName,".mat"),"wb");
  if(fp==NULL)
	 printf("File could not be opened.\n");
  
  else{
	 x.type = 0000;				  /* 0000  in pc:s (and in alphas? ) */
	 x.mrows = rows;
	 x.ncols = cols;
	 x.imagf = isImag;;
	 x.namelen = strlen(matrixName)+1;
	 fwrite(&x,sizeof(Fmatrix),1,fp);
	 fwrite(matrixName, sizeof(char), x.namelen,fp);
	 mn = x.mrows *x.ncols;
	 fwrite(realData,sizeof(double),mn,fp);
	 if(x.imagf)
		fwrite(imagData,sizeof(double),mn,fp);
  }
  
  fclose(fp);
  return 0;
}
