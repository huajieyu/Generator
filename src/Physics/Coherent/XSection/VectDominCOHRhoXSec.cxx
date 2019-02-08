//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - March 11, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 24, 2016 - Libo Jiang, Chris Kullenburgh, Sanjib Mishra
   Was renamed to VectDominCOHRhoXSec
For the rho-nucleon total cross section,

Return the coherent rho cross section in GeV^-2
calculated from the formula given by
Kopliovich and Marage: Int. Journal of Modern Physics A
Vol 8, No 9(1993) 1512-1602 Eq 156

Rhosig is the cross section of the pi-N scattering
given by Equation 2.16 in Coherent Diffractive Pion Production
in Charged Current Neutrino Interactions. Lyle John Winton.
and Compilation of cross section I-pi and Pi+ induced reactions
Technical Report 79-01, CERN-HERA, 1979


 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/VectDominCOHRhoXSec.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/GHEP/GHepParticle.h"
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
VectDominCOHRhoXSec::VectDominCOHRhoXSec() :
XSecAlgorithmI("genie::VectDominCOHRhoXSec")
{

}
//____________________________________________________________________________
VectDominCOHRhoXSec::VectDominCOHRhoXSec(string config) :
XSecAlgorithmI("genie::VectDominCOHRhoXSec", config)
{

}
//____________________________________________________________________________
VectDominCOHRhoXSec::~VectDominCOHRhoXSec()
{

}
//____________________________________________________________________________
double VectDominCOHRhoXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();


  double RhoMass = interaction->ExclTag().RhoMass();
  bool is_NC  = interaction->ProcInfo().IsWeakNC();  
  
  //----- Compute the coherent rho^0 production d2xsec/dxdy
  //      see page 34 in Nucl.Phys.B223:29-144 (1983)
  double E      = init_state.ProbeE(kRfLab); // neutrino energy
  double x      = kinematics.x(); // bjorken x
  double y      = kinematics.y(); // inelasticity y
  double Q2     = kinematics.Q2();  // momentum transfer Q2>0
  double w      = kinematics.W();
  double w2     = w*w;
  // Libo added for coherent rho
  //double v      = k - 0.5 * q2/Mnuc;
  //double v2     = TMath::Power(v, 2);
  double ELEP   = (1.0-y)*E; // rho's energy
  //calculate the variables needed for coherent rho, Libo
  double ENU2   = E*E;
  double RFAC   = 0;
  double RFLG1  = 1;
  double RFLG2  = 0;
  double UU     = E*y;
  double UU2    = UU*UU;
  double Q      = TMath::Sqrt(UU*UU+Q2);  
  //initialize variables for the test of cross-section
  /*E=1;
  UU=1;
  UU2=1;
  ELEP=1;
  y=1;
  Q=1;
  Q2=1;
  */
  //--------------------------------------------------
  
  double FABS;
  //temp declare
  double AMT = interaction->InitState().Tgt().Mass();   //nucleus mass
  double AMT2=AMT*AMT;
  double ATM=(double) init_state.Tgt().A();    //get the atm number
  double ATM2=ATM*ATM;   //

  Q2=2.0*AMT*E*x*y;
  Q=TMath::Sqrt(UU*UU+Q2);

  w2=2.0*AMT*UU-Q2+AMT2;
  w=TMath::Sqrt(w2);

  double RA=R0*TMath::Power(ATM, 1.0/3.0); 
  double BSLOPE=(1.0/3.0)*TMath::Power((RA*FMTOGEV),2.0);
  //calculate of TMIN and TMAX

  //evrec->FindParticle(213,kIStStableFinalState,0); 
  
  //GHepParticle *p=evrec->FindParticle(213,kIStStableFinalState,0);

  //double AMMES=p->Mass();

  double AMMES=RhoMass;
  
  double AMMES2 = AMMES*AMMES;
  //std::cout<<"[VectDomin] "<<RhoMass<<std::endl;

  double EN2C=(w2+AMT2-AMMES2)/(2.0*w);
  double PN2C=TMath::Sqrt(EN2C*EN2C-AMT2);
  double EN1C=(w2+AMT2+Q2)/(2*w);
  double PN1C=TMath::Sqrt(EN1C*EN1C-AMT2);

  double TMIN=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C-PN1C)*(PN2C-PN1C);
  double TMAX=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C+PN1C)*(PN2C+PN1C);
  //------------------------------------------------
  double S0=(AMT+AMMES)*(AMT+AMMES);
  //std::cout<<"[VectDomin] Ev= "<<E<<" UU= "<<UU<<" Elep = "<<ELEP<<" Q2= "<<Q2<<" TMIN= "<<TMIN<<" TMAX= "<<TMAX<<std::endl;

  //std::cout<<"S0= "<<S0<<std::endl; 
  //Initialize some of the variables for the test of cross-section
  //TMIN=1;
  //TMAX=2;
  //----------------------------------------------------  
  double BTNUM   = exp(-BSLOPE*TMIN)-exp(-BSLOPE*TMAX);
  double EPSILON=(4.0*E*ELEP-Q2)/(4.0*E*ELEP+Q2+2.0*UU2);
 
  double PROPRHO=AMRHO2/(Q2+AMRHO2);


  double GRHO2=2.4*4.0*kPi;
  double RHOSIG=0.0;
  double ICOHPIN=2;
  if(ICOHPIN==1){
    RHOSIG=24.0;
  }
  else if(ICOHPIN==2){
  RHOSIG=24.0*(1.0+0.5/sqrt(UU));
  }
  
  double R=RFAC*(RFLG1+RFLG2*Q2/AMRHO2);
  
  // calculate the rho obsorption factor in the nucleus
  if(R>1.0) {R=1.0;}
  if(ATM<1.1){
  FABS=1.0;
  }
  else{
  FABS=exp(-9.0/16.0/kPi*RA/(pow(R0, 3.0))/FM2TOMB*19.0*(1.0+0.44891/sqrt(UU)));
  }

  //for the test of the cross-section calculation
  

  //------------------------------------------
  double SIGOUT = GFERMI2/(2.0*kPi*kPi)*
        1.0/GRHO2*
        PROPRHO*PROPRHO*
        Q/ENU2*
        Q2/(1.0 - EPSILON)*
        (1.0 + EPSILON*R)*
        ATM2/(16.0*kPi)*
        (RHOSIG*RHOSIG)*
        FABS*
        BTNUM/BSLOPE*
        MBTOCM2/GEVTOMB;
  SIGOUT=SIGOUT/GEVTOMB/MBTOCM2; //convert the unit from 10^-40cm^2 to GeV^-2
  //multiply by Jacobian factor to translate from dQ2dv to dxdy
  SIGOUT=2.0*AMT*y*ENU2*SIGOUT;
  if(isnan(SIGOUT)){SIGOUT=0;}

  //if(is_NC) {SIGOUT=SIGOUT*0.5*(1.0-2.0*sinthetaw2)*(1.0-2.0*sinthetaw2);}
  //std::cout<<"[VECTDOMIN]"<<"is NC???"<<is_NC<<SIGOUT<<std::endl;
  double xsec=SIGOUT;
  //===========================================================
  // effect of pion absorption in the nucleus

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("VectDominCohRho", pDEBUG)
      << "\n momentum transfer .............. Q2    = " << Q2
      << "\n mass number .................... A     = " << A
      << "\n Rho mass .................... AMMES   = " << AMMES
      << "\n propagator term ................ proprho = " << PROPRHO
      << "\n Re/Im of fwd pion scat. ampl. .. Re/Im = " << fReIm
      << "\n total pi+N cross section ....... sigT  = " << sTot
      << "\n inelastic pi+N cross section ... sigI  = " << sInel
      << "\n nuclear size scale ............. R0    = " << R0
      << "\n rho absorption factor ......... Fabs  = " << FABS
      << "\n t integration range ............ [" << TMIN << "," << TMAX << "]"
      << "\n t integration factor ........... tint  = " << tint;
#endif

  // compute the cross section for the CC case

  /*  //only consider charged current here
  if(interaction->ProcInfo().IsWeakCC()) { 
     // Check whether a modification to Adler's PCAC theorem is applied for
     // including the effect of the muon mass. 
     // See Rein and Sehgal, PCAC and the Deficit of Forward Muons in pi+ 
     // Production by Neutrinos, hep-ph/0606185
     double C = 1.;
     if(fModPCAC) {
        double ml    = interaction->FSPrimLepton()->Mass();
        double ml2   = TMath::Power(ml,2);
        double Q2min = ml2 * y/(1-y);
        if(Q2>Q2min) {
           double C1    = TMath::Power(1-0.5*Q2min/(Q2+kPionMass2), 2);
           double C2    = 0.25*y*Q2min*(Q2-Q2min)/ TMath::Power(Q2+kPionMass2,2);
           C = C1+C2;
        } else {
           C = 0.;
        }
     }
     xsec *= (2.*C); 
  }*/


#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalCohRho", pINFO)
         << "d2xsec/dxdy[COHPi] (x= " << x << ", y="
                       << y << ", E=" << E << ") = "<< xsec;
#endif

  //----- The algorithm computes d^2xsec/dxdy
  //      Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double VectDominCOHRhoXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool VectDominCOHRhoXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  int nu = init_state.ProbePdg();

  if (!proc_info.IsCoherent())  return false;
  if (!proc_info.IsWeak())      return false;
  if (target.HitNucIsSet())     return false;
  if (!(target.A()>1))            return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void VectDominCOHRhoXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VectDominCOHRhoXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VectDominCOHRhoXSec::LoadConfig(void)
{
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();

//talk with Gabe about the parameters for coherent rho
/*  fMa      = fConfig->GetDoubleDef("Ma",            gc->GetDouble("COH-Ma"));
  fReIm    = fConfig->GetDoubleDef("Re-Im-Ampl",    gc->GetDouble("COH-ReImAmpl"));
  fRo      = fConfig->GetDoubleDef("Ro",            gc->GetDouble("COH-Ro"));
  fModPCAC = fConfig->GetBoolDef("UseModifiedPCAC", gc->GetBool("COH-UseModifiedPCAC"));
*/
  GetParam( "COH-Ma",fMa ) ;
  GetParam( "COH-Ro", fRo ) ;

  
  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

