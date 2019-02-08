//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Renamed COHPiXSec -> COHRhoXSec. Adapt to naming changes made to the coherent 
   generator for including coherent vector meson production.
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/XSection/COHRhoXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/Utils/Genfun.h"
//#include "Framework/Utils/cohrhosigxy.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace genie::controls;
//____________________________________________________________________________
COHRhoXSec::COHRhoXSec() :
XSecIntegratorI("genie::COHRhoXSec")
{

}
//____________________________________________________________________________
COHRhoXSec::COHRhoXSec(string config) :
XSecIntegratorI("genie::COHRhoXSec", config)
{

}
//____________________________________________________________________________
COHRhoXSec::~COHRhoXSec()
{

}
//____________________________________________________________________________
double COHRhoXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHRhoXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  double AMT = in->InitState().Tgt().Mass();   //nucleus mass 
  double AMT2=AMT*AMT;
  //double ATM=(double) init_state.Tgt().A();
  double AMMES=0.7700;
  //genie::utils::genfun::xgencal(2*kPionMass, 2*AMRHO,  &AMMES);  

  in->ExclTagPtr()->SetRhoMass(AMMES);
  double AMMES2=AMMES*AMMES;  
  //std::cout<<"[COHRHOXSEC]"<<in->ExclTag().RhoMass()<<std::endl;
  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  //const EventGeneratorI * evg = rtinfo->RunningThread();
  //const XSecAlgorithmI *model = evg->CrossSectionAlg();
  
  
  int ibadev=0;
  double xmax, ymin, ymax, gx_temp, gy_temp;
  double xsec=0.0;
  double total_xsec=0.0;
  double W2, W, Q2, Elep, UU;
  int nmax=1000000;


  RandomGen * rnd = RandomGen::Instance();
  if(Ev<1.0) {return 0;}
  else{

  for(int k=0; k<nmax; k++){

  int Nthtry=0;
  do{
  gy_temp = rnd->RndKine().Rndm();
  gx_temp = rnd->RndKine().Rndm();


  Nthtry=Nthtry+1;
  ibadev=0;

  double Amlep=kMuonMass;
  double Amlep2=Amlep*Amlep;

  double S0=(AMT+AMMES)*(AMT+AMMES);
  xmax=(AMT2+2.0*AMT*Ev*gy_temp-S0)/(2.0*AMT*Ev*gy_temp);

  if(gx_temp>xmax) {ibadev=1;}
  Elep=Ev*(1.0-gy_temp);
  if(Elep<Amlep){ibadev=2;}

  double Pklep=TMath::Sqrt(Elep*Elep-Amlep2);

  double AA=1.0-Pklep/Elep;
  double BB=1.0+Pklep/Elep;
  ymax=(Ev*BB-Amlep2/(2.0*Ev))/(Ev*BB+AMT*gx_temp);
  ymin=(Ev*AA-Amlep2/(2.0*Ev))/(Ev*AA+AMT*gx_temp);

  if(gy_temp>ymax){ibadev=3;}
  if(gy_temp<ymin){ibadev=4;}

  UU=Ev*gy_temp;
  Q2=2.0*AMT*UU*gx_temp;
  double Q=TMath::Sqrt(UU*UU+Q2);

  W2=2.0*AMT*UU-Q2+AMT2;
  W=TMath::Sqrt(W2);

  if(W2<=S0){ibadev=5;}
  if(Pklep<=0. || Q<=0.){
  ibadev=6;
  }
  } while(ibadev !=0); //end of do loop


  double EN2C=(W2+AMT2-AMMES2)/(2.0*W);
  double PN2C=sqrt(EN2C*EN2C-AMT2);

  double EN1C=(W2+AMT2+Q2)/(2.0*W);
  double PN1C=sqrt(EN1C*EN1C-AMT2);

  //cout<<"WW2="<<WW2<<"AMT2="<<AMT2<<"AMMES2="<<AMMES2<<"W="<<W<<endl;
  //cout<<"EN1C= "<<EN1C<<" "<<"EN2C= "<<EN2C<<endl;
  //cout<<"PN1C= "<<PN1C<<" "<<"PN2C= "<<PN2C<<endl;
  //calculate TMIN and TMAX Here for test
  double TMIN=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C-PN1C)*(PN2C-PN1C);
  double TMAX=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C+PN1C)*(PN2C+PN1C);

  //TMIN2=((AMMES2+QQ2)/(2.0*UU))*((AMMES2+QQ2)/(2.0*UU));


  interaction->KinePtr()->Setx(gx_temp);
  interaction->KinePtr()->Sety(gy_temp);
  interaction->KinePtr()->SetQ2(Q2);
  interaction->KinePtr()->SetW(W); 
  
  xsec=0.0;
  if(gx_temp>0 && gy_temp>0){
   //std::cout<<"[COHRhoXSec] Ev= "<<Ev<<" UU= "<<UU<<" Elep = "<<Elep<<" Q2= "<<Q2<<" TMIN= "<<TMIN<<" TMAX= "<<TMAX<<std::endl;
   if (model->Id().Name() == "genie::VectDominCOHRhoXSec") {

   //-- Access cross section algorithm for running thread
   //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
   //const EventGeneratorI * evg = rtinfo->RunningThread();
   //const XSecAlgorithmI *model = evg->CrossSectionAlg();

   xsec = model->XSec(interaction, kPSxyfE);
   //xsec=genie::utils::cohrhosigxy::cohrhosigxy_cal(Ev, UU, Elep, Q2, TMIN, TMAX);   
   total_xsec +=xsec;
   } //end of if xsec model
  } //end of gx_temp>0 and gy_temp>0  
  } //end of loop over k


  xsec=total_xsec/(nmax*1.0);
  } // end of loop over else energy>1.0
  //const InitialState & init_state = in->InitState();
  //double Ev = init_state.ProbeE(kRfLab);
  LOG("COHRhoXSec", pINFO) << "XSec[COH] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void COHRhoXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHRhoXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHRhoXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Get GSL integration type & relative tolerance

  GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1E-2 ) ;

  int max_eval, min_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  GetParamDef( "gsl-min-eval", min_eval, 5000 ) ;

  fGSLMaxEval  = (unsigned int) max_eval ;
  fGSLMinEval  = (unsigned int) min_eval ;


  //-- COH model parameter t_max for t = (q - p_pi)^2
  //fTMax = fConfig->GetDoubleDef("COH-t-max", gc->GetDouble("COH-t-max"));
  //-- COH model bounds of integration for Q^2
  //fQ2Min = fConfig->GetDoubleDef("COH-Q2-min", gc->GetDouble("COH-Q2-min"));
  //fQ2Max = fConfig->GetDoubleDef("COH-Q2-max", gc->GetDouble("COH-Q2-max"));
}
//____________________________________________________________________________
