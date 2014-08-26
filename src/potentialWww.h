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


/*! 
  \brief initializes module at startup
  Initializes the module. par and pos have to be initialized at this point.
*/
void initwwwPot(struct parameters *par,struct systemPos pos);

/*! 
  \brief Calculates energy of system. 
  Calculates the total energy of the whole system pos.enPot is updated
  \param calcForces If true then also the forces in pos.xa,pos.ya,pos.za are updated
  
 */
double wwwPot(struct parameters *par,struct systemPos pos,int calcForces);

double partwwwPot(struct parameters *par,struct systemPos pos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold);
void reinitWwwPotAfterStep(struct parameters *par,struct systemPos pos,int stepAccepted);
void reinitWwwAfterBondChange(struct parameters *par,struct systemPos pos);
