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




static double keatingInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  int  j,k;
  double rik2,rij2;
  
  int numOfNgbrs;
  int listStart,listEnd;
  int jlpos,klpos;
 
  double forcej,forcek;
  double dotProd;
  double cosTerm;
  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;
  
  struct bondList *blist;
  getBondList(&blist);
  
  numOfNgbrs=numBonds(i);
  listStart=blist->head[i];
  listEnd=listStart+numOfNgbrs;
  
  for(jlpos=listStart;jlpos<listEnd;jlpos++){
    j=blist->list[jlpos]; /*j is the ngbr*/
  
    drx_ij=pos.x[j]-pos.x[i]; 
    dry_ij=pos.y[j]-pos.y[i]; 
    drz_ij=pos.z[j]-pos.z[i];	
    rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
	 
    if(rij2>par->minhBox2){
      PERIODIC(drx_ij,dry_ij,drz_ij);
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
    }
	 
    /*first the pair part*/
    /*here we sum over  i=0...N-1, j>i*/
      
    if( j>i){ /*dont doublecount!! in pair part*/
      poti=KSTRETCH(i,j)*POW2(R02(i,j)-rij2);
      totPot+=poti;
      pos.enPot[i]+=0.5*poti;
      pos.enPot[j]+=0.5*poti;
		
      if(calcForces){ /*j=n*/
		  forcej=4.0*KSTRETCH(i,j)*(R02(i,j)-rij2);
		  /*xyza is acceleration, NOT force*/
		  pos.xa[i] -= forcej*drx_ij;
		  pos.ya[i] -= forcej*dry_ij;
		  pos.za[i] -= forcej*drz_ij;
			  
		  pos.xa[j] += forcej*drx_ij;
		  pos.ya[j] += forcej*dry_ij;
		  pos.za[j] += forcej*drz_ij;
      }
    }
		 
		 
    /*then the three body part */
    /*here we sum over i=0...N-1, j!=i, k!=i && k>j */
    for(klpos=jlpos+1;klpos<listEnd;klpos++){		 
      k=blist->list[klpos]; /*k is the ngbr*/		
     
      drx_ik=pos.x[k]-pos.x[i]; 
      dry_ik=pos.y[k]-pos.y[i]; 
      drz_ik=pos.z[k]-pos.z[i];	
		
      rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
      if(rik2>par->minhBox2){
		  PERIODIC(drx_ik,dry_ik,drz_ik);
		  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
      }
		
      dotProd=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij);
      cosTerm=(dotProd-COS0(i)*R0(i,j)*R0(i,k));
      poti=KBEND(j,i,k)*POW2(cosTerm);
      pos.enPot[i]+=poti;
      totPot+=poti;

		
      if(calcForces){
		  cosTerm*=-2.0*KBEND(j,i,k);

		  forcej=cosTerm*drx_ik;
		  forcek=cosTerm*drx_ij;
		  pos.xa[j] += forcej;
		  pos.xa[k] += forcek;
		  pos.xa[i] += -forcej-forcek;
		  
		  forcej=cosTerm*dry_ik;
		  forcek=cosTerm*dry_ij;	 
		  pos.ya[j] += forcej;
		  pos.ya[k] += forcek;
		  pos.ya[i] += -forcej-forcek;
			
		  forcej=cosTerm*drz_ik;
		  forcek=cosTerm*drz_ij;  
		  pos.za[j] += forcej;
		  pos.za[k] += forcek;
		  pos.za[i] += -forcej-forcek;
      }
    }
  }
  return totPot;
}



static double simpKeatingInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  
  double rij,rij2;
  double rik,rik2;
  
  int numOfNgbrs;
  int listStart,listEnd;
  int  j,k;
  int jlpos,klpos;
  
  double forcej,forcek;
  double varTemp1,varTemp2,varTemp3,varTemp4;
  double cosijk,cosTerm;
  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;

  struct bondList *blist;

  getBondList(&blist);
  
  numOfNgbrs=numBonds(i);
  listStart=blist->head[i];
  listEnd=listStart+numOfNgbrs;
    
  for(jlpos=listStart;jlpos<listEnd;jlpos++){
    j=blist->list[jlpos]; /*j is the ngbr*/
    
    drx_ij=pos.x[j]-pos.x[i]; 
    dry_ij=pos.y[j]-pos.y[i]; 
    drz_ij=pos.z[j]-pos.z[i];	
    rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
	 
    if(rij2>par->minhBox2){
      PERIODIC(drx_ij,dry_ij,drz_ij);
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
    }
    rij=sqrt(rij2);
    /*first the pair part*/
    /*here we sum over  i=0...N-1, j>i*/
	 
    if( j>i){ /*dont doublecount!! in pair part*/
      poti=KSTRETCH(i,j)*POW2(R0(i,j)-rij);
      totPot+=poti;
      pos.enPot[i]+=0.5*poti;
      pos.enPot[j]+=0.5*poti;
	  
      if(calcForces){ /*j=n*/
		  forcej=2.0*KSTRETCH(i,j)*(R0(i,j)-rij)/rij;
		  /*xyza is acceleration, NOT force*/
		  pos.xa[i] -= forcej*drx_ij;
		  pos.ya[i] -= forcej*dry_ij;
		  pos.za[i] -= forcej*drz_ij;
			  
		  pos.xa[j] += forcej*drx_ij;
		  pos.ya[j] += forcej*dry_ij;
		  pos.za[j] += forcej*drz_ij;
      }
    }
	 
		 
    /*then the three body part */
    /*here we sum over i=0...N-1, j!=i, k!=i && k>j */
    for(klpos=jlpos+1;klpos<listEnd;klpos++){		 
      k=blist->list[klpos]; /*k is the ngbr*/
      drx_ik=pos.x[k]-pos.x[i]; 
      dry_ik=pos.y[k]-pos.y[i]; 
      drz_ik=pos.z[k]-pos.z[i];	
		
      rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
      if(rik2>par->minhBox2){
		  PERIODIC(drx_ik,dry_ik,drz_ik);
		  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
      }
      rik=sqrt(rik2);
		
      cosijk=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/(rik*rij);
      cosTerm=(cosijk-COS0(i));
		
		poti=KBEND(j,i,k)*POW2(cosTerm);

      pos.enPot[i]+=poti;
      
      totPot+=poti;
		
		  
      if(calcForces){
		
		  varTemp1=-2.0*KBEND(j,i,k)*cosTerm;

		  varTemp2=varTemp1/(rik*rij);
		  varTemp3=varTemp1*cosijk/rij2;
		  varTemp4=varTemp1*cosijk/rik2;
		  
		  
		  forcej=varTemp2*drx_ik-varTemp3*drx_ij;
		  forcek=varTemp2*drx_ij-varTemp4*drx_ik;
		  pos.xa[j] += forcej;
		  pos.xa[k] += forcek;
		  pos.xa[i] += -forcej-forcek;
		  
		  forcej=varTemp2*dry_ik-varTemp3*dry_ij;
		  forcek=varTemp2*dry_ij-varTemp4*dry_ik;
		  pos.ya[j] += forcej;
		  pos.ya[k] += forcek;
		  pos.ya[i] += -forcej-forcek;
		  
		  forcej=varTemp2*drz_ik-varTemp3*drz_ij;
		  forcek=varTemp2*drz_ij-varTemp4*drz_ik;
		  pos.za[j] += forcej;
		  pos.za[k] += forcek;
		  pos.za[i] += -forcej-forcek;
      }
    }
  }
  return totPot;
}

/*remember energy units not taken care of!*/

static double stixInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces){
  double drx_ij,dry_ij,drz_ij;
  double drx_ik,dry_ik,drz_ik;
  double drx_jk,dry_jk,drz_jk;
  int  j,k;
  double rik2,rij2,rjk2,rjk,rik,rij;
  
  int numOfNgbrs;
  int listStart,listEnd;
  int jlpos,klpos;

  double forcei,forcej,forcek;
  double dotProd;
  double cosjik,dcosijk;
  double alphajik,dalphajik;

  double totPot=0.0;						  /* here energy is calculated from beginning */
  double poti;
  double fcij,fcik,dfcij,dfcik;
  double fcmax;
  double expbetaij,expgammaij;
  double expbetaik,expgammaik;
  
  double tempTerm1,tempTerm2,tempTerm3,tempTerm4;
  
  struct bondList *blist;
  getBondList(&blist);

  
  expgammaij=exp(wwwPar.stix_gamma*(1.62-wwwPar.stix_rc));
  fcij=1.0/(1.0+expgammaij);                           
  dfcij=0.0;
  expgammaik=exp(wwwPar.stix_gamma*(1.62-wwwPar.stix_rc));
  fcik=1.0/(1.0+expgammaik);                           
  dfcik=0.0;
  
  numOfNgbrs=numBonds(i);
  listStart=blist->head[i];
  listEnd=listStart+numOfNgbrs;
  
  /*start potential for Si atoms (first two terms)*/
  if( pos.aType[i]==SI_ATYPE)
	 for(jlpos=listStart;jlpos<listEnd;jlpos++){
		j=blist->list[jlpos]; /*j is the ngbr*/
		
		drx_ij=pos.x[j]-pos.x[i]; 
		dry_ij=pos.y[j]-pos.y[i]; 
		drz_ij=pos.z[j]-pos.z[i];	
		rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		if(rij2>par->minhBox2){
		  PERIODIC(drx_ij,dry_ij,drz_ij);
		  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		}
		rij=sqrt(rij2);
		/*first part of pot*/
		
		/*
		  expgammaij=exp(wwwPar.stix_gamma*(rij-wwwPar.stix_rc));
		  fcij=1.0/(1.0+expgammaij);                          
		*/
		expbetaij=exp(wwwPar.stix_beta*(rij-wwwPar.stix_r0));
		poti=fcij*wwwPar.stix_D*expbetaij*(expbetaij-2);
		totPot+=poti;
      pos.enPot[i]+=0.5*poti;
      pos.enPot[j]+=0.5*poti;

		if(calcForces){ 
		  /*		  dfcij=-fcij*fcij*expgammaij*wwwPar.stix_gamma;*/
		  forcei=expbetaij*wwwPar.stix_D/rij*(dfcij*(expbetaij-2)+fcij*2*wwwPar.stix_beta*(expbetaij-1));
	
		  pos.xa[i] += forcei*drx_ij;
		  pos.ya[i] += forcei*dry_ij;
		  pos.za[i] += forcei*drz_ij;
			  
		  pos.xa[j] -= forcei*drx_ij;
		  pos.ya[j] -= forcei*dry_ij;
		  pos.za[j] -= forcei*drz_ij;
      }
		
		
		/*then the three body part */
		/*here we sum over i=0...N-1, j!=i, k!=i && k>j */
		for(klpos=jlpos+1;klpos<listEnd;klpos++){		 
		  k=blist->list[klpos]; 
		  drx_ik=pos.x[k]-pos.x[i]; 
		  dry_ik=pos.y[k]-pos.y[i]; 
		  drz_ik=pos.z[k]-pos.z[i];	
		  
		  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
		  if(rik2>par->minhBox2){
			 PERIODIC(drx_ik,dry_ik,drz_ik);
			 rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
		  }
		  rik=sqrt(rik2);
		  cosjik=(drx_ik*drx_ij+dry_ik*dry_ij+drz_ik*drz_ij)/(rik*rij);
		  alphajik=acos(cosjik);
		  
		  
		  if(rik>rij){ /*if rik <rij we do not need to calculate fcik*/
			 /*
				expgammaik=exp(wwwPar.stix_gamma*(rik-wwwPar.stix_rc));
				fcik=1.0/(1.0+expgammaik);                           
				
			 */
			 fcmax=fcik;
		  }
		  else
			 fcmax=fcij;
		  
		  poti=fcmax*wwwPar.stix_galpha*POW2(alphajik-wwwPar.stix_alpha0);

		  
		  totPot+=poti;
		  pos.enPot[i]+=poti;
		  if(calcForces){
			 
			 tempTerm1=fcmax*wwwPar.stix_galpha*2*(alphajik-wwwPar.stix_alpha0)/sqrt(1-POW2(cosjik))/(rij*rik);
			 tempTerm2=tempTerm1*(-cosjik);
			 
			 forcek=tempTerm1*drx_ij+tempTerm2*drx_ik*rij/rik;
			 forcej=tempTerm1*drx_ik+tempTerm2*drx_ij*rik/rij;
			 pos.xa[k]+=forcek;
			 pos.xa[j]+=forcej;
			 pos.xa[i]+=-forcek-forcej; 
			 forcek=tempTerm1*dry_ij+tempTerm2*dry_ik*rij/rik;
			 forcej=tempTerm1*dry_ik+tempTerm2*dry_ij*rik/rij;
			 pos.ya[k]+=forcek;
			 pos.ya[j]+=forcej;
			 pos.ya[i]+=-forcek-forcej; 
			 forcek=tempTerm1*drz_ij+tempTerm2*drz_ik*rij/rik;
			 forcej=tempTerm1*drz_ik+tempTerm2*drz_ij*rik/rij;
			 pos.za[k]+=forcek;
			 pos.za[j]+=forcej;
			 pos.za[i]+=-forcek-forcej; 
			 
			 if(rij>rik){
				
				tempTerm4=-dfcij*wwwPar.stix_galpha*POW2(alphajik-wwwPar.stix_alpha0)/rij;
				forcej=tempTerm4*drx_ij;
				pos.xa[j]+=forcej;
				pos.xa[i]+=-forcej; 
				forcej=tempTerm4*dry_ij;
				pos.ya[j]+=forcej;
				pos.ya[i]+=-forcej; 	
				forcej=tempTerm4*drz_ij;
				pos.za[j]+=forcej;
				pos.za[i]+=-forcej; 
			 }
			 else{
				/*				dfcik=-fcik*fcik*expgammaik*wwwPar.stix_gamma;*/
				tempTerm4=-dfcik*wwwPar.stix_galpha*POW2(alphajik-wwwPar.stix_alpha0)/rik;
				forcek=tempTerm4*drx_ik;
				pos.xa[k]+= forcek;
				pos.xa[i]+=-forcek; 
				forcek=tempTerm4*dry_ik;
				pos.ya[k]+= forcek;
				pos.ya[i]+=-forcek; 	
				forcek=tempTerm4*drz_ik;
				pos.za[k]+= forcek;
				pos.za[i]+=-forcek; 
				
			 }
		  }
		  
		  
		}
	 }
  /*start poetntial for Oxygen atoms*/
  if(pos.aType[i]==O_ATYPE)
	 for(jlpos=listStart;jlpos<listEnd;jlpos++){
		j=blist->list[jlpos]; /*j is the ngbr*/
		
		drx_ij=pos.x[j]-pos.x[i]; 
		dry_ij=pos.y[j]-pos.y[i]; 
		drz_ij=pos.z[j]-pos.z[i];	
		rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		if(rij2>par->minhBox2){
		  PERIODIC(drx_ij,dry_ij,drz_ij);
		  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		}
		rij=sqrt(rij2);

		/*first part of pot*/
		/*
		  expgammaij=exp(wwwPar.stix_gamma*(rij-wwwPar.stix_rc));
		  fcij=1.0/(1.0+expgammaij);                           
		*/
		/*then the three body part */
		/*here we sum over i=0...N-1, j!=i, k!=i && k>j */
		
		for(klpos=jlpos+1;klpos<listEnd;klpos++){		 
		  k=blist->list[klpos]; 
		  drx_ik=pos.x[k]-pos.x[i]; 
		  dry_ik=pos.y[k]-pos.y[i]; 
		  drz_ik=pos.z[k]-pos.z[i];	
		  
		  rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
		  if(rik2>par->minhBox2){
			 PERIODIC(drx_ik,dry_ik,drz_ik);
			 rik2=LENGTH2(drx_ik,dry_ik,drz_ik);
		  }
		  rik=sqrt(rik2);
		  
		  drx_jk=drx_ij-drx_ik;
		  dry_jk=dry_ij-dry_ik;
		  drz_jk=drz_ij-drz_ik;
		  rjk2=LENGTH2(drx_jk,dry_jk,drz_jk);
		  rjk=sqrt(rjk2);

		  if(rik>rij){ /*if rik <rij we do not need to calculate fcik*/
			 /*			 
							expgammaik=exp(wwwPar.stix_gamma*(rik-wwwPar.stix_rc));
							fcik=1.0/(1.0+expgammaik);                           
			 */
			 fcmax=fcik;
		  }
		  else
			 fcmax=fcij;
		  		  
		  poti=fcmax*wwwPar.stix_gL*POW2(rjk-wwwPar.stix_L0);
		  
		  totPot+=poti;
		  pos.enPot[j]+=0.5*poti;
		  pos.enPot[k]+=0.5*poti;

		  if(calcForces){
			 tempTerm1=fcmax*wwwPar.stix_gL*2*(rjk-wwwPar.stix_L0)/rjk;
			 pos.xa[k]+=tempTerm1*drx_jk;
			 pos.xa[j]-=tempTerm1*drx_jk;
			 pos.ya[k]+=tempTerm1*dry_jk;
			 pos.ya[j]-=tempTerm1*dry_jk;
			 pos.za[k]+=tempTerm1*drz_jk;
			 pos.za[j]-=tempTerm1*drz_jk;

			 if(rij>rik){
				/*
				  dfcij=-fcij*fcij*expgammaij*wwwPar.stix_gamma;
				*/
				tempTerm1=-dfcij*wwwPar.stix_gL*POW2(rjk-wwwPar.stix_L0)/rij;
				pos.xa[j]+=tempTerm1*drx_ij;
				pos.xa[i]-=tempTerm1*drx_ij; 
				pos.ya[j]+=tempTerm1*dry_ij;
				pos.ya[i]-=tempTerm1*dry_ij; 
				pos.za[j]+=tempTerm1*drz_ij;
				pos.za[i]-=tempTerm1*drz_ij; 
			 }
			 else{
				/*				dfcik=-fcik*fcik*expgammaik*wwwPar.stix_gamma;*/
				tempTerm1=-dfcik*wwwPar.stix_gL*POW2(rjk-wwwPar.stix_L0)/rik;
				pos.xa[k]+=tempTerm1*drx_ik;
				pos.xa[i]-=tempTerm1*drx_ik; 
				pos.ya[k]+=tempTerm1*dry_ik;
				pos.ya[i]-=tempTerm1*dry_ik; 
				pos.za[k]+=tempTerm1*drz_ik;
				pos.za[i]-=tempTerm1*drz_ik; 
				
			 }
		  }
		}
	 }
  
  /*and now the last non-bonded repulsive part*/
  if(  pos.aType[i]==O_ATYPE){
	 
	 int ngbrListPosj;
 
	 
	 ngbrListPosj=nd.head[i]; /*where the ngbrs start*/
	 while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
		j=nd.list[ngbrListPosj++]; /*j is the ngbr*/
		
		if(pos.aType[j]!=O_ATYPE) 
		  continue;

		drx_ij=pos.x[j]-pos.x[i]; 
		dry_ij=pos.y[j]-pos.y[i]; 
		drz_ij=pos.z[j]-pos.z[i];	
		rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		if(rij2>par->minhBox2){
		  PERIODIC(drx_ij,dry_ij,drz_ij);
		  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
		}
		rij=sqrt(rij2);
		/*
		  poti=wwwPar.stix_A*exp(-rij/wwwPar.stix_b);
		*/

		poti= tabulatedVal(rij,wwwPar.stix_repTable);
		/*	printf("%g %g %g\n",rij,poti,wwwPar.stix_A*exp(-rij/wwwPar.stix_b));*/
		totPot+=poti;
      pos.enPot[i]+=0.5*poti;
      pos.enPot[j]+=0.5*poti;
		
		
		
		if(calcForces){ 
		  /*forcei=-poti/(wwwPar.stix_b*rij);*/
		  forcei=tabulatedDeriv(rij,wwwPar.stix_repTable)/rij;
		  
		  pos.xa[i] += forcei*drx_ij;
		  pos.ya[i] += forcei*dry_ij;
		  pos.za[i] += forcei*drz_ij;
		  pos.xa[j] -= forcei*drx_ij;
		  pos.ya[j] -= forcei*dry_ij;
		  pos.za[j] -= forcei*drz_ij;
      }
		
		
	 }
  
	 
  }
  
  return totPot;
}





static double tabulatedRepulsiveInnerLoop(int i,struct parameters *par,struct systemPos pos,int calcForces){
  int  j,k;
  int jlpos,klpos;
  
  double drx,dry,drz,r2,r;
  
  double totPot,poti;
  double force;
  struct bondList *blist;
  int numOfNgbrs;
  int listStart,listEnd;
  struct valueTable *tempTable;
  
  

	 
  totPot=0.0;  
  getBondList(&blist);
  
	
  numOfNgbrs=numBonds(i);
  listStart=blist->head[i];
  listEnd=listStart+numOfNgbrs;
  
  for(jlpos=listStart;jlpos<listEnd;jlpos++){
	 j=blist->list[jlpos]; /*j is the ngbr*/
	 for(klpos=jlpos+1;klpos<listEnd;klpos++){		 
		k=blist->list[klpos]; /*k is the ngbr*/
		/*choose pottable type*/
		if(pos.aType[j]==O_ATYPE && pos.aType[k]==O_ATYPE)
		  tempTable=&wwwPar.ooPottable;
		else  if(pos.aType[j]==SI_ATYPE && pos.aType[k]==SI_ATYPE)
		  tempTable=&wwwPar.ssPottable;
		else
		  tempTable=&wwwPar.soPottable;
		
		drx=pos.x[k]-pos.x[j]; 
		dry=pos.y[k]-pos.y[j]; 
		drz=pos.z[k]-pos.z[j];	
		
		r2=LENGTH2(drx,dry,drz);
		if(r2>par->minhBox2){
		  PERIODIC(drx,dry,drz);
		  r2=LENGTH2(drx,dry,drz);
		}
		r=sqrt(r2);
		poti= tabulatedVal(r,*tempTable);
		pos.enPot[j]+=0.5*poti;
		pos.enPot[k]+=0.5*poti;
		totPot+=poti;
		if(calcForces){
		  force=-tabulatedDeriv(r,*tempTable);
		  pos.xa[j] -= force*drx/r;
		  pos.xa[k] += force*drx/r;
			  
		  pos.ya[j] -= force*dry/r;
		  pos.ya[k] += force*dry/r;
		  
		  pos.za[j] -= force*drz/r;
		  pos.za[k] += force*drz/r;
		}
		
	 }
  }
  
  
  return totPot;
}



static double repulsiveInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces){
  int j;
  int ngbrListPosj;
  double poti,totPot;
  double forcej;
  double drx_ij,dry_ij,drz_ij;
  double rij2;
  	 
  totPot=0.0;
  ngbrListPosj=nd.head[i]; /*where the ngbrs start*/
  while(nd.list[ngbrListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
    j=nd.list[ngbrListPosj++]; /*j is the ngbr*/
    
    if(i>j){ /*dont doublecount!! in pair part*/
      drx_ij=pos.x[j]-pos.x[i]; 
      dry_ij=pos.y[j]-pos.y[i]; 
      drz_ij=pos.z[j]-pos.z[i];	
      rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      if(rij2>par->minhBox2){
		  PERIODIC(drx_ij,dry_ij,drz_ij);
		  rij2=LENGTH2(drx_ij,dry_ij,drz_ij);
      }
      if(rij2<wwwPar.repR02){
		  poti=wwwPar.repK*POW2(wwwPar.repR02-rij2);
		  totPot+=poti;
		  pos.enPot[i]+=0.5*poti;
		  pos.enPot[j]+=0.5*poti;

		  if(calcForces){ /*j=n*/
			 forcej= 4.0*wwwPar.repK*(wwwPar.repR02-rij2);
			 
			 /*xyza is acceleration, NOT force*/
			 pos.xa[i] -= forcej*drx_ij;
			 pos.ya[i] -= forcej*dry_ij;
			 pos.za[i] -= forcej*drz_ij;
			 
			 pos.xa[j] += forcej*drx_ij;
			 pos.ya[j] += forcej*dry_ij;
			 pos.za[j] += forcej*drz_ij;
		  }	
      }
    }
    
  }

  return totPot;

}

double suboxideInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces){
  double poti=0.0;
  int numberOfOngbrs,bondListPosj;
  int j;
  struct bondList *blist;
  getBondList(&blist);
  
  if(pos.aType[i]==siType){
    bondListPosj=blist->head[i]; /*where the ngbrs start*/
    numberOfOngbrs=0; 
    while(blist->list[bondListPosj]!=-1){ /*get ngbrs of i, stops when the ngbrs of atom i+1 comes*/
      j=blist->list[bondListPosj]; /*j is the ngbr*/
      bondListPosj++;
      if(pos.aType[j]==oType)
		  numberOfOngbrs++;		  
    }
    poti=wwwPar.subOxidePenalty[numberOfOngbrs];
    pos.enPot[i]+=poti;
  }
  
  
  return poti;
}
 

static double danglBondsInnerLoop(int  i,struct parameters *par,struct systemPos pos,int calcForces){
  double poti=0.0;
  struct bondList *blist;
  getBondList(&blist);

  if(pos.aType[i]==siType)
    poti=blist->danglBonds[i]*1.0*EV;
  else if(pos.aType[i]==oType)
    poti=blist->danglBonds[i]*4.0*EV;
	 
  pos.enPot[i]+=poti;
  return poti;
} 
