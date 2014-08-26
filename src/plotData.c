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
#include <string.h>
#include "plotData.h"

static char fname[10];
/*one could later on put more variables so the user of the module could have more control of the plots*/

void initPlotData(plotData* pd,char label[]){
  static int num=0;
 
  
  /*lets create the commands*/
  num++;
  sprintf(fname,"%s%d%s","data",num,".dat");
  
  pd->plotCommand=malloc(256*sizeof(char));
  pd->programStart=malloc(256*sizeof(char));
  pd->setValues=malloc(256*sizeof(char));
  
  
  sprintf(pd->setValues,"%s\"%s\"\n","set terminal x11;set nokey;set title",label);
  sprintf(pd->plotCommand,"%s%s%s\n","plot '",fname,"' with lines\n");
  strcpy(pd->programStart,"/usr/bin/gnuplot>tempdump\n");
  
  pd->plotFile=fopen(fname,"w");
  if( pd->plotFile==NULL){
	 printf("Error can't temp file for plotdata\n");
	 exit(2);
  }
  
  pd->plotPipe=popen(pd->programStart, "w"); /*open pipe to gnuplot*/
  
  if (pd->plotPipe == NULL){
	 printf("could not find gnuplot");
	 exit(2);
  }
  
  fprintf(pd->plotPipe, "%s", pd->setValues); /*init terminal*/ 
  fflush(pd->plotPipe);

}
void clearPlotData(plotData* pd){
  fclose(pd->plotFile);
  pd->plotFile=fopen(fname,"w");
}
void destroyPlotData(plotData* pd){
  free(pd->plotCommand);
  free(pd->programStart);
  free(pd->setValues);
  pclose(pd->plotPipe);
  fclose(pd->plotFile);
  /*fclose...*/
} 
void putPlotVal(plotData* pd,double x,double y){
  fprintf(pd->plotFile,"%g %g\n",x,y);
  fflush(pd->plotFile); 
}
void plotNow(plotData* pd){
  fprintf(pd->plotPipe, "%s", pd->plotCommand);
  fflush(pd->plotPipe);
} 
