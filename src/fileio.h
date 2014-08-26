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


#ifndef FILEIO
#define FILEIO

struct fileIo{
  char inFileName[256];   /*the name of the file from where we read in the system*/
  char inFileFormat[256]; /*its type*/
  
  char inBondFile[256];   /*the name of the file from where we read in the bonds of the system*/
  
  char outFileName[256];  /*the name of the file where we write in the system*/
  char outFileFormat[256];/*its type*/
  
  char comment[1024];     /*a comment of the system*/

};
void initStructFileIo(struct fileIo *fio);




void openSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par);
void readBondFile(struct fileIo *fio,struct systemPos *pos,struct parameters *par);

void writeSystem(struct fileIo *fio,struct systemPos *pos,struct parameters *par, char comment[]);
void writeSystemOpendx(struct fileIo *fio,struct systemPos *pos,struct parameters *par);

void openParameterFile(char *filename);
void closeParameterFile(void);
int getParValue(char *parameter,void *value,char *formatString);


#endif
