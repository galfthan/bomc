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


#include "ngbrs.h"


void mcVolume(struct systemPos pos,struct parameters *par,int restore);
void backupState(struct systemPos pos); /* this saves some backup information about current state, after some change has been made and energy minimized this can be used to save this state.   */
void restoreState(struct systemPos pos); /* restore previously backuped state (backups position, energy  and bonds which have been changed using one of the mc steps*/

double accProbNPT(struct parameters *par,double potVal,double oldPotVal,int nAtoms,double oldVol);
double accProbNVT(struct parameters *par,double potVal,double oldPotVal);

int mcBondDiffuse(struct systemPos *pos,int *part,int *atomsInPart);
int mcBondSwitch(struct systemPos *pos ,int *part,int *atomsInPart);
int makeSpecificSwitch(struct systemPos *pos,int ai,int aj, int aib, int ajb); /* just makes the switch if it is at all possible, does not look at how close they are orwhat kind of rings are born */
int mcBondCreate(struct systemPos *pos ,struct parameters *par,  int *part,int *atomsInPart);
int mcBondBreak(struct systemPos *pos  ,int *part,int *atomsInPart);
int mcOxygenDiffuse(struct systemPos *pos, int *part,int *atomsInPart);


double getKT(int iteration,int nAtoms);
