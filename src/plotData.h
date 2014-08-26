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


#ifndef PLOTDATA_DEFS
#define PLOTDATA_DEFS

#define GNUPLOT "/usr/bin/gnuplot"


typedef struct plotData
{
  FILE *plotPipe; 
  FILE *plotFile;
  char *plotCommand ;
  char *programStart ;
  char *setValues ;
} plotData;

void initPlotData(plotData* pd,char label[]);
void destroyPlotData(plotData* pd);
void plotNow(plotData* pd);
void putPlotVal(plotData* pd,double x,double y);
void clearPlotData(plotData* pd);

  
#endif

