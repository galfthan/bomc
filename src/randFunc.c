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
#define LEN 100
#define MAXLONG 2147483648.0


static int numOfSeries;
static int *seed;  
static double *randNums;

void cr250(int *locSeed,double *locRandNums);


void initRandNum(int locNumOfSeries,long givenSeed){
  static int firsttime=1;
  int i;
  
  if(firsttime){
    numOfSeries=locNumOfSeries;
    seed=malloc(sizeof(int)*numOfSeries*(LEN+250));
    randNums=malloc(sizeof(double)*numOfSeries*LEN);
  }


  srandom(givenSeed);
  
  for(i=0;i<numOfSeries*(LEN+250);i++)
    seed[i]=random(); 
  firsttime=0;
}

double randNum(int series ){
  static int firsttime=1;
  static int *currRandNum;
  
  if(firsttime){
	 int i;
	 currRandNum=malloc(sizeof(int)*numOfSeries);
	 for(i=0;i<numOfSeries;i++){
		cr250(seed+i*(LEN+250),randNums+i*LEN);
		currRandNum[i]=0;
	 }
	 firsttime=0;
  }
  
  if(currRandNum[series]==LEN){
    cr250(seed+series*(LEN+250),randNums+series*LEN);
    currRandNum[series]=0;
  }
  return (randNums+series*LEN)[currRandNum[series]++];
}






/* R250 Pseudorandom number generator c-version by Tuukka Sales
   (sales@rock.helsinki.fi) 24.4.1998 - tested to give the same
   sequence as the F77 version.

   This program generates uniform pseudorandom numbers between (0,1]

   */

void cr250(int *locSeed,double *locRandNums)
{
  int i;
  int q=103, p=250;
  for(i=0; i<LEN; i++){
    locSeed[i+p]=locSeed[i+p-q]^locSeed[i];
    locRandNums[i]=(double)(locSeed[i+p]/MAXLONG);
  }
  for(i=0; i<p; i++)
    locSeed[i]=locSeed[i+LEN];
}






