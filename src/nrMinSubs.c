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

#include <stdlib.h>
#include <stdio.h>
#include "shared.h"
#include "nrMinSubs.h"
#include "initialize.h"
#include <memory.h>
#include <math.h>




#define TOL 0.1
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ITMAX 50
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

/*Numerical recipies subroutines*/
/*changes made to linmin and f1dim, mnbrak and brent left alone..*/

/*these globals are needed for f1dim*/
static struct systemPos comPos;
static struct systemPos *comOrgPos;
static struct parameters *comPar;
double *comDirx,*comDiry,*comDirz;

int numOfPotCalcs;
static double (*comPartPot)(struct parameters*,struct systemPos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold);
static double (*comPot)(struct parameters*,struct systemPos,int);
static int *comPart;
static int  comAtomsInPart;

double brent(double ax, double bx, double cx,double fbx,double (*f)(double), double tol, double *xmin);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,double *fc, double (*func)(double));
double f1dim(double x);
double f1dimPart(double x);

void linmin(struct systemPos *pos,double *dir,struct parameters *par,double *fret,int allowAdaptive, double (*pot)(struct parameters*,struct systemPos,int))
{
  static int firsttime=1;
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  static double meanXmin=-6.0e-6;

  /*these can point to the same place*/
  comPar=par; 
  comOrgPos=pos;
  comDirx=dir;
  comDiry=dir+pos->nAtoms;
  comDirz=dir+2*pos->nAtoms;
  comPot=pot;
  /*compos has to be a copy, it is changed in f1dim*/
  if(firsttime){
    initStructPos(&comPos,pos->nAtoms);
    memcpy(comPos.aType,pos->aType,sizeof(int)*pos->nAtoms);
    memcpy(comPos.aFlag,pos->aFlag,sizeof(long)*pos->nAtoms);
    firsttime=0;
  }
  
  numOfPotCalcs=0;
  /*the minimization in itself*/
  ax=meanXmin/2.0; //-1.0e-6;
  xx=meanXmin;     //-6.0e-6
  fx=f1dim(xx);
  if(allowAdaptive && fx<(*fret)){  /* lower energy, accept this step directly*/
    xmin=xx;
    *fret=fx; 
  }
  else{
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=brent(ax,xx,bx,fx,f1dim,TOL,&xmin);
    meanXmin=0.8*meanXmin+0.2*xmin; /* adaptive, 1% comes from the calculated xmin*/
  }
  /*update system to the next step*/
  for (j=0;j<pos->nAtoms;j++) {
    pos->x[j]+=comDirx[j]*xmin;
    pos->y[j]+=comDiry[j]*xmin;
    pos->z[j]+=comDirz[j]*xmin;
  }
  
}

void linminPart(struct systemPos *pos,double *dir,struct parameters *par,double *fret,int allowAdaptive, int *part, int atomsInPart
		,double (*partpot)(struct parameters*,struct systemPos,int *part, int *atomsInPart,int calcForces,int updatePart,int returnRealPot,double updateTreshold))
{
  static int firsttime=1;
  int j,jps;
  double xx,xmin,fx,fb,fa,bx,ax;
  static double meanXmin=-6.0e-6;

  /*these can point to the same place*/
  comPar=par; 
  comOrgPos=pos;
  comDirx=dir;
  comDiry=dir+pos->nAtoms;
  comDirz=dir+2*pos->nAtoms;
  comPartPot=partpot;
  comPart=part;
  comAtomsInPart=atomsInPart;

  /*compos has to be a copy, it is changed in f1dim*/
  if(firsttime){
    initStructPos(&comPos,pos->nAtoms);
    memcpy(comPos.aType,pos->aType,sizeof(int)*pos->nAtoms);
    memcpy(comPos.aFlag,pos->aFlag,sizeof(long)*pos->nAtoms);
    
    firsttime=0;
  }
  
  /*positions cant be assumed to be identical.. each time this is run, f1doim will update only part of them */
  
  memcpy(comPos.x,pos->x,sizeof(double)*pos->nAtoms);
  memcpy(comPos.y,pos->y,sizeof(double)*pos->nAtoms);
  memcpy(comPos.z,pos->z,sizeof(double)*pos->nAtoms);

  
  numOfPotCalcs=0;
  /*the minimization in itself*/
  ax=meanXmin/2.0; //-1.0e-6;
  xx=meanXmin;     //-6.0e-6
  fx=f1dimPart(xx);
  if(allowAdaptive && fx<(*fret)){		 	    
    xmin=xx;
    *fret=fx; 
  }
  else{
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dimPart);
    *fret=brent(ax,xx,bx,fx,f1dimPart,TOL,&xmin);
    meanXmin=0.99*meanXmin+0.01*xmin; /* adaptive, 1% comes from the calculated xmin*/
  }
  // printf("xx %g\n",meanXmin);
 
  /*update system to the next CG step*/
  for(jps=0;jps<atomsInPart;jps++){
    j=part[jps];
    pos->x[j]+=comDirx[j]*xmin;
    pos->y[j]+=comDiry[j]*xmin;
    pos->z[j]+=comDirz[j]*xmin;
  
  }
}


void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double)){
  double ulim,u,r,q,fu,dum;
  
  *fa=(*func)(*ax);
  //*fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,(*func)(u))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}


double brent(double ax, double bx, double cx,double fbx, double (*f)(double), double tol,
	     double *xmin)
{
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=fbx;
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }

  printf("Too many iterations in brent\n");
  
  return fx; /*lets go in some direction at least.. */
}


double f1dim(double x){
  int j;

  numOfPotCalcs++;
  /*update system to the next step*/
  for (j=0;j<comPos.nAtoms;j++){
    comPos.x[j]=comOrgPos->x[j]+comDirx[j]*x;
    comPos.y[j]=comOrgPos->y[j]+comDiry[j]*x;
    comPos.z[j]=comOrgPos->z[j]+comDirz[j]*x;
   
  }

  return (*comPot)(comPar,comPos,0);
}


double f1dimPart(double x){
  int j,jps;

  numOfPotCalcs++;
  /*update system to the next step*/
 
  for(jps=0;jps<comAtomsInPart;jps++){
    j=comPart[jps];
    comPos.x[j]=comOrgPos->x[j]+comDirx[j]*x;
    comPos.y[j]=comOrgPos->y[j]+comDiry[j]*x;
    comPos.z[j]=comOrgPos->z[j]+comDirz[j]*x;
   
  }
  
  return (*comPartPot)(comPar,comPos,comPart,&comAtomsInPart,0,0,0,0.0);
}
