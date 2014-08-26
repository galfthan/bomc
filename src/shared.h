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


#ifndef SHARED
#define SHARED


/*==== MACROS ====*/

/*the PERIODIC macros assume we have  par-> */

/*#define MAKEPERIODIC(r,hlength,length) MAKEPERIODIC_DEF_DIM(r,-hlength,hlength,length);*/
/*#define MAKEPERIODIC(r,hlength,length) ( -rint(r/length)*length)*/
/*#define PERIODICXYZ(r,num) if(par->periodic[num]) MAKEPERIODIC(r,par->hBox[num],par->box[num]);*/


#define MAKEPERIODIC_DEF_DIM(r,left,right,delta) (r>right)?r-=delta:( (r<left)?r+=delta:0.0 );
#define MAKEPERIODIC(r,length) (r-=rint(r/length)*length)
//#define SYSTEMISPERIODIC 


//#define PERIODICXYZ(r,num) if(par->periodic[num]) { (r+par->hBox[num])/ }
#define PERIODICXYZ(r,num) if(par->periodic[num]) MAKEPERIODIC_DEF_DIM(r,-par->hBox[num],par->hBox[num],par->box[num]);
/*#define PERIODICXYZ(r,num) if(par->periodic[num]) MAKEPERIODIC(r,par->box[num]);*/

#define PERIODIC(rx,ry,rz) PERIODICXYZ(rx,0); PERIODICXYZ(ry,1); PERIODICXYZ(rz,2);

//#define PERIODIC(rx,ry,rz) PERIODICXYZ(rx,0); PERIODICXYZ(ry,1); PERIODICXYZ(rz,2);

#define POW2(x) ((x) * (x))
#define POW3(x) (POW2(x)*(x))
#define POW4(x) (POW2(x)*POW2(x))
#define POW5(x) (POW3(x)*POW2(x))
#define POW6(x) (POW3(x)*POW3(x))
#define POW7(x) (POW3(x)*POW4(x))
#define POW8(x) (POW4(x)*POW4(x))
#define POW9(x) (POW5(x)*POW4(x))
#define POW10(x) (POW5(x)*POW5(x))


#define MAX(x,y) ((x)<(y)?(y):(x))
#define MAX3(x,y,z) (MAX((x),(y))>(z)?MAX((x),(y)):(z))

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MIN3(x,y,z) (MIN((x),(y))<(z)?MIN((x),(y)):(z))

#define LENGTH2(x,y,z) (POW2(x)+POW2(y)+POW2(z))

#define VOLUME (par->box[0]*par->box[1]*par->box[2])

#define DEBUG 0


//#define FLAG(f,i) (( (f) >> (i))&1)
//#define ATOM_IN_SI(f) FLAG(f,0)


/*======UNITS========*/

/*=1ps, in seconds */
#define UNITTIME 1.0e-12 

/*=1Å in meters  */  
#define UNITDIST 1.0e-10 

/*=1u, in kg*/
#define UNITMASS 1.66053873e-27   

//#define UNITCHARGE 4.074972637944948e-17 /*in coulomb, =sqrt(mass*length)*length/time  */
/*with this unitcharge the energy in the coulomb interactions will be in the correct units*/
//#define ETOUNITCHARGE 3.93175005663e-3 /*ECHARGE/UNITCHARGE*/

/*=1e, in coulomb*/
#define ECHARGE 1.60217739e-19 
#define UNITCHARGE ECHARGE
#define ETOUNITCHARGE 1.0


#define TOEV 1.036426266132741e-04
/*= mass*(length/time)^2/echarge */
/*this is the value of the energy in EV when the energy is 1 in the UNITS defined above*/
/*based on the UNITS defined above, make sure it is correct!!!*/

#define EV 9.648539724213475e+03
/*TOEV^-1.0 */
/*this is the value EV in the UNITS  defined above, */
/*based on the UNITS defined above, make sure it is correct!!!*/
           
#define PASCAL 6.022141982800968e-08/*=U_time^2*U_dist/U_mass*/
#define TOPASCAL 16.605387e6

#define JOULETOEV  6.2415097e+18
#define EVTOJOULE  ECHARGE


#define TOJOULE 1.66053873e-23
#define JOULE 6.022141982800966e+22


#define KCALTOEV 	2.6131953e+22


/*========CONSTANTS ========*/

 
#define EPSILON0_SI 8.8541878e-12 
/*Vs/Am*/
#define KB_JOULE 	1.3806503e-23 
/*joule/k */
#define KBEV 8.617342239757910e-05 
/*in ev/K */
#define AVOGADRO 6.022142e23

#define AMORPH_FLAG 0
#define FIXED_FLAG 1
#define CRYST_FLAG 2
#define INTERFACE_FLAG 3



#define SI_ATYPE 1
#define O_ATYPE 0


#define FALSE 0
#define TRUE 1
#define NUMOFATYPES 13



/*num of possible atomic types, increase if neccessary*/

/*======= STRUCTS =======*/
struct systemPos
{
  int    nAtoms;         /*number of atoms in system*/
  double *x,*y,*z;       /*position of atoms*/  
  double *prevx,*prevy,*prevz; /*the positions of the atoms after the last accepted step*/
  double *xa,*ya,*za;    /*force on each atom*/
  int    *aType;         /*number which tells of which material the atom is
									atypes are:
									0  1  
									O  Si
								 */
  double *enPot; /*potential  energy of each atom*/
  double *prevEnPot; /*contains the potentialEnergy of the system after the last accepted step (corresponding to prevxyz*/
  long *aFlag;           /*addtitional information about each atom, see 
									FLAG macro*/
};

typedef struct partStruct{
  int n;
  int *list;
  int *table;
} sysPart;


struct parameters{
  int periodic[3];                    /*tells if the system is periodic in x,y,z direction */
  int varVol[3];
  double box[3],hBox[3],minhBox2; /*the size and half the size of the system in A and the min of the hBoxes squared*/
  double volume;
  double skin;                        /*skin (used by potentials) in A*/
  char aName[NUMOFATYPES][4];         /*tells the name of a certain atype*/
  double aMass[NUMOFATYPES];          /*tells the mass of a certain atype (in internal units)*/ 
  double kTrand;
  double kT;
  double pressure;
  long *seed;                         /*seed number for random number generator*/
  int interfaceSys;				  /* tells if we have a system with a interface between si/siO2 . if 1 then bond switches are not performed in Si area and simAnneal is biased towards siO2 */
  
};



#endif
