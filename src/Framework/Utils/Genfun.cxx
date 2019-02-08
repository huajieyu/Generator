//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 22, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Framework/Utils/Genfun.h"
//#include "Conventions/Constants.h"

#include "fgrhom.h"

using namespace genie;
//using namespace genie::constants;

//______________________________________________________________________
double genie::utils::genfun::xgencal(
                       double XMI, double XMA,  double* XGEN)
{
//----------------------------------------------------------------------
int INIT, NT, II; 
double X, DX, XL, FL, FGEN;
double FNT[101], FTOT;
NT=100;
INIT=0;

//cout<<"XMI= "<<XMI<<" "<<"XMA= "<<XMA<<endl;
INIT=1;
DX=(XMA-XMI)/double (NT);
X=XMI+DX/2.0;
FNT[1]=fgrhom(X)*DX;

//get the total area of the fgrhom
for(int I=2; I<=NT; I++){
   X=X+DX;
   FNT[I]=FNT[I-1]+fgrhom(X)*DX;
}
//integral of the bw function
FTOT=FNT[NT];

//generate a random number between 0 and FTOT;
FGEN=((double) rand() / double (RAND_MAX))*FTOT;

int I=0;
do {
I=I+1;
II=I;
}while(FNT[I]<FGEN);

     XL=XMI+(II-1)*DX;
     FL=0.;
     if(II>=2) FL=FNT[II-1];
     *XGEN=XL+DX*(FGEN-FL)/(FNT[II]-FL);
  

return 0;



//-----------------------------------------------------------------------




}
//______________________________________________________________________

