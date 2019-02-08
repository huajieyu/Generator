//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 Authors: Libo Jiang, Sanjib Mishra, Chris Kullenburgh
 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (EVGModules 
   package)
 @ Mar 03, 2009 - CA
   Renamed COHPiKinematicsGenerator -> COHRhoKinematicsGenerator in
   anticipation of reusing the code for simulating coherent production of
   vector mesons.
 @ May 06, 2009 - CA
   Fix a problem with the search for the max cross section over the allowed
   phase space which prevented kinematics to be generated for events near the 
   energy threshold.
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <TMath.h>
#include <TFile.h>
#include <TF2.h>
#include <TNtupleD.h>
#include <TSystem.h>
#include <TROOT.h>
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/EventGen/COHRhoKinematicsGenerator.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Genfun.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHRhoKinematicsGenerator::COHRhoKinematicsGenerator() :
KineGeneratorWithCache("genie::COHRhoKinematicsGenerator")
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHRhoKinematicsGenerator::COHRhoKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::COHRhoKinematicsGenerator", config)
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHRhoKinematicsGenerator::~COHRhoKinematicsGenerator()
{
  if(fEnvelope) delete fEnvelope;
}
//___________________________________________________________________________
void COHRhoKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("COHRhoKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant

  //-------------------------------------------------------------------
  //get the mass of the rho meason
  double AMMES;
  genie::utils::genfun::xgencal(2*kPionMass, 2*AMRHO,  &AMMES);  
  //AMMES=0.7700; //set a constant for test
  interaction->ExclTagPtr()->SetRhoMass(AMMES);
  
  int  probepdgc = interaction->InitState().ProbePdg();
  if(probepdgc==14){ 
  interaction->ExclTagPtr()->SetRhoPDG(kPdgRhoP);
  }else if(probepdgc==-14){
  interaction->ExclTagPtr()->SetRhoPDG(kPdgRhoM);
  }

  /*
        
  */    

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1.;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
      this->throwOnTooManyIterations(iter,evrec);
     }


     const InitialState & init_state = interaction->InitState();
     double Ev    = init_state.ProbeE(kRfLab);
     double AMT = interaction->InitState().Tgt().Mass();   //nucleus mass
     //double AMT=11.17793; //set a constant for test
     double AMT2=AMT*AMT;
     double ATM=(double) init_state.Tgt().A();    //get the atm number
 
     //---------------------------------------------------------------------
     // call the function ComputeXYLim here to get the XY Limits and gx gy 
     Range1D_t varx, vary;
     double gx, gy, gt, gQ2, gW, gt_temp;
     unsigned int tryxy=0;
     int index_xy=0;
     label5: 
     COHRhoKinematicsGenerator::ComputeXYLim(interaction, AMMES, &varx, &vary, &gx, &gy, &gt, &gQ2, &gW,  &index_xy); 
     tryxy=index_xy;
     if(tryxy>=kRjMaxIterations)
     {this->throwOnTooManyIterations(iter,evrec); }

     //---------------------------------------------------------------------
     // Reset the kps.YLim() and kps.XLim() here
     const KPhaseSpace & kps = interaction->PhaseSpace();

     Range1D_t y = kps.YLim();
     Range1D_t x = kps.XLim();
 
     x.min=varx.min;
     x.max=varx.max;
     y.min=vary.min;
     y.max=vary.max;
 

     //--------------------------------------------------------------------

     double v=Ev*gy;

     
     //------------------------------------------------------------------
     //-- decide whether to accept the current kinematics
    
     interaction->KinePtr()->Sett(gt);
     interaction->KinePtr()->SetQ2(gQ2);
     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);  
     interaction->KinePtr()->SetW(gW);
     std::cout<<"[COHRhoKinematicsGenerator] hadronic invariance mass is "<<gW<<std::endl;
     // computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);
            
     const Kinematics &   kinematics = interaction -> Kine();

     //-- For the subsequent kinematic selection with the rejection method:
     //   Calculate the max differential cross section or retrieve it from the
     //   cache. Throw an exception and quit the evg thread if a non-positive
     //   value is found.
     //   If the kinematics are generated uniformly over the allowed phase
     //   space the max xsec is irrelevant


    double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
    //-------------------------------------------------------------------  
    
    int NCOHNUC=4;
    int COHLOCH=4;
    int INUClibo=0;
    int IElibo=0;
    double COHXSEC[301][5];
    double COHXMAX[301][7];
    const int NBINE=300;

   
    /*
    if(gSystem->Getenv("GENIE")){
    std::string cohrho0=string(gSystem->Getenv("GENIE"));
    std::string cohrho1=cohrho0+string("/src/Physics/Coherent/EventGen/cohxsec.dat");
    if(!(gSystem->AccessPathName(cohrho1.c_str()))){
      LOG("COH", pINFO)<<"Load COHRHO data from: "<<cohrho1;
      std::ifstream file(cohrho1.c_str());
      std::vector<int> numbers;
      std::string line;
      
      for(int INUC=1; INUC<=NCOHNUC; INUC++){
         for(int IE=1; IE<=NBINE; IE++){ 
         std::getline(file, line);
         std::string raw_nums;
         std::istringstream iss(line);
         //while(std::getline(iss, raw_nums, ' ')){
         //   numbers.push_back(std::strod(raw_nums.c_str()));
         //}
         //Spline * spline = new Spline(INUClibo,IElibo,COHXSEC[IE][INUC],COHXMAX[IE][INUC]);
         //std::stringstream s_str0(line0);
         iss>>INUClibo>>IElibo>>COHXSEC[IE][INUC]>>COHXMAX[IE][INUC];
         if(INUC==COHLOCH){
           COHXMAX[IE][NCOHNUC+1]=COHXMAX[IE][INUC];
           COHXMAX[IE][NCOHNUC+2]=COHXMAX[IE][INUC];
         }
         //cohrho1>>INUClibo>>IElibo>>COHXSEC[IE][INUC]>>COHXMAX[IE][INUC];
         }
      }      
    }
    }
    
    

    
    
    const double  DE=1.0;
    int newIE=(Ev/DE)+1;
    int ICOHNORM=1; 
    double TOPTEST=0.0;
    double XSECTEST=0.0;
    if(newIE<=NBINE) {
       TOPTEST=COHXMAX[newIE][ICOHNORM];
       XSECTEST=COHXSEC[newIE][ICOHNORM];
    }

    else{
        TOPTEST=COHXMAX[NBINE][ICOHNORM];
        XSECTEST=COHXSEC[NBINE][ICOHNORM];
    }
    */
     
     accept = true; 
     
     //std::cout<<TOPTEST/GEVTOMB/MBTOCM2<<std::endl;
     
     if(xsec_max * rnd->RndKine().Rndm() > xsec) {goto label5;}

     //-------------------------------------------------------------------
 
  
     double Q2     = kinematics.Q2();  // momentum transfer Q2>0
     double W      = kinematics.W(); 
     AMMES= interaction->ExclTag().RhoMass();
     double AMMES2=AMMES*AMMES; 
     double W2=W*W;
     double EN2C=(W2+AMT2-AMMES2)/(2.0*W);
     double PN2C=sqrt(EN2C*EN2C-AMT2);

     double EN1C=(W2+AMT2+Q2)/(2.0*W);
     double PN1C=sqrt(EN1C*EN1C-AMT2);

     double TMIN=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C-PN1C)*(PN2C-PN1C);
     double TMAX=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C+PN1C)*(PN2C+PN1C);

     double RA=R0*TMath::Power(ATM, 1.0/3.0);
     double BSLOPE=(1.0/3.0)*TMath::Power((RA*FMTOGEV),2.0);
     do{
     gt_temp = - (1.0/BSLOPE)*log((double) rand() / double (RAND_MAX));
     }while(gt_temp<TMIN ||gt_temp>TMAX);


     interaction->KinePtr()->Sett(gt_temp);
    

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("COHRhoKinematics", pNOTICE) << "Selected: x = "<< gx << ", y = "<< gy;

        // the COH cross section should be a triple differential cross section
        // d^2xsec/dxdydt where t is the the square of the 4p transfer to the
        // nucleus. The cross section used for kinematical selection should have
        // the t-dependence integrated out. The t-dependence is of the form
        // ~exp(-bt). Now that the x,y kinematical variables have been selected
        // we can generate a t using the t-dependence as a PDF.
        //get the value of AMMES
     
 
        //const Kinematics &   kinematics = interaction -> Kine();


        Q2=2.0*AMT*v*gx;
        double Q=TMath::Sqrt(v*v+Q2);
        
        W2=2.0*AMT*v-Q2+AMT2;
        W=TMath::Sqrt(W2);
        AMMES2=AMMES*AMMES;     

        EN2C=(W2+AMT2-AMMES2)/(2.0*W);
        PN2C=sqrt(EN2C*EN2C-AMT2);

        EN1C=(W2+AMT2+Q2)/(2.0*W);
        PN1C=sqrt(EN1C*EN1C-AMT2);

        TMIN=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C-PN1C)*(PN2C-PN1C); 
        TMAX=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C+PN1C)*(PN2C+PN1C); 

        RA=R0*TMath::Power(ATM, 1.0/3.0);
        BSLOPE=(1.0/3.0)*TMath::Power((RA*FMTOGEV),2.0);
 
        gt = kinematics.t();  
        double PNf[4], PMES[4],PKMES, PKNf;
        double PLEP[4], cosal, sinal,PKLEP, costh, sinth, PTLEP, PLLEP, ELEP;
        double cosx1,sinx1;
        //generate gt between TMIN and TMAX from a random number
        bool nonacceptgt=true;
       
        unsigned int ind=0;
        do {
         ind++;
         if(ind > kRjMaxIterations) {
          this->throwOnTooManyIterations(ind,evrec);
         }
         PLEP[0]=Ev-v;
         ELEP=PLEP[0];
         PKLEP=TMath::Sqrt(PLEP[0]*PLEP[0]-kMuonMass2);
         costh=(2*Ev*PLEP[0]-kMuonMass2-Q2)/(2*Ev*PKLEP);
         if(TMath::Abs(costh)>1.0) {costh=1.0; sinth=0.0;}
         sinth=TMath::Sqrt(1.0-costh*costh);

         PNf[0]=AMT+gt/(2.0*AMT);
         PKNf=TMath::Sqrt(PNf[0]*PNf[0]-AMT2);

         PMES[0]=v+AMT-PNf[0];
         PKMES=TMath::Sqrt(PMES[0]*PMES[0]-AMMES2);

         cosx1=(Q2-gt-AMMES2+2*PMES[0]*v)/(2.0*Q*PKMES);
         if(TMath::Abs(cosx1)>1.0) {cosx1=1.0; sinx1=0.0;}
         sinx1=TMath::Sqrt(1.0-cosx1*cosx1);
  
         PLLEP=PKLEP*costh;
         PTLEP=PKLEP*sinth;
         sinal=-PTLEP/Q;
         if(TMath::Abs(sinal)>1.0) {sinal=1.0; cosal=0.0;}
         cosal=TMath::Sqrt(1.0-sinal*sinal);
 
          if (TMath::Abs(cosx1)<1.0 && TMath::Abs(costh)<1.0 && TMath::Abs(sinal)<1.0){
           nonacceptgt=false;
          }else{nonacceptgt=true;}
        } while (nonacceptgt==true);


        double b=BSLOPE;
        double tsum=(TMath::Exp(-b*TMIN) - TMath::Exp(-b*TMAX))/b;

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->Setx(gx, true);
        interaction->KinePtr()->Sety(gy, true);
        interaction->KinePtr()->Sett(gt_temp, true);
        interaction->KinePtr()->SetW(W, true);
        interaction->KinePtr()->SetQ2(Q2, true);
        interaction->KinePtr()->ClearRunningValues();

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec*TMath::Exp(-b*gt)/tsum,kPSxytfE);
        /*
        double LAM=2.0*kPi*((double) rand() / double (RAND_MAX));
        cosla=cos(LAM);
        sinla=sin(LAM);
        PLEP[1]=PTLEP*cosla;
        PLEP[2]=PTLEP*sinla;
        PLEP[3]=PLLEP;
        PLEP[0]=PLEP[0];


        // get the rho meson's 4 momentum
        PLMES=PKMES*cosx1;
        PTMES=PKMES*sinx1;
        double PHI=2.0*kPi*((double) rand() / double (RAND_MAX));
        sinph=sin(PHI); cosph=cos(PHI);
        PMES[1]=PTMES*sinph;
        PMES[2]=PTMES*cosph;
        PMES[3]=PLMES;

        Q2=2.0*AMT*v*gx;
        Q=TMath::Sqrt(v*v+Q2);



        double PLNUC=Q-PLMES;
        double PTNUC=-PTMES;

        PNf[1]=PTNUC*sinph;
        PNf[2]=PTNUC*cosph;
        PNf[3]=PLNUC; 
        
        //COHRhoKinematicsGenerator::dgdrot(PMES, cosal, sinal, cosla, sinla);
        //COHRhoKinematicsGenerator::dgdrot(PNf,  cosal, sinal, cosla, sinla);
       
        double P1=PMES[1];
        double P2=PMES[2];
        double P3=PMES[3];
        PMES[1]=P1*cosal*cosla-P2*sinla+P3*sinal*cosla;
        PMES[2]=P1*cosal*sinla+P2*cosla+P3*sinal*sinla;
        PMES[3]=-P1*sinal+P3*cosal;

        double PP1=PNf[1];
        double PP2=PNf[2];
        double PP3=PNf[3];
        PNf[1]=PP1*cosal*cosla-PP2*sinla+PP3*sinal*cosla;
        PNf[2]=PP1*cosal*sinla+PP2*cosla+PP3*sinal*sinla;
        PNf[3]=-PP1*sinal+PP3*cosal;
        */      

        //-------------------------------------------------------

       //-------------------------------------------------------
        return;
     }
  }// iterations
}
//___________________________________________________________________________
double COHRhoKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHRhoKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif
  double max_xsec = 0.;
  double AMMES = in->ExclTag().RhoMass();  
  const Kinematics &   kinematics = in -> Kine();
  const InitialState & init_state = in -> InitState();


  double Ev = init_state.ProbeE(kRfLab);

  const int Nx = 100;
  const int Ny = 100;
  //need a tolerance factor? not sure
  //double Q2 = kinematics.Q2();
  double ggy  = kinematics.y();
  double ggx  = kinematics.x();
  double AMT = in->InitState().Tgt().Mass();   //nucleus mass
  double AMT2=AMT*AMT;
  double S0 = (AMT+AMMES)*(AMT+AMMES);  
  Range1D_t x, y;
  x.max=(AMT2+2.0*AMT*Ev*ggy-S0)/(2.0*AMT*Ev*ggy);
  double Elep=Ev*(1.0-ggy);
  double Amlep=kMuonMass;
  double Amlep2=Amlep*Amlep;
  double Pklep=TMath::Sqrt(Elep*Elep-Amlep2);
  double AA=1.0-Pklep/Elep;
  double BB=1.0+Pklep/Elep;
  y.max =(Ev*BB-Amlep2/(2.0*Ev))/(Ev*BB+AMT*ggx);
  y.min =(Ev*AA-Amlep2/(2.0*Ev))/(Ev*AA+AMT*ggx);
  x.min = kASmallNum;

  //--------------------------------------------------
  //const KPhaseSpace & kps = in->PhaseSpace();
  //Range1D_t y = kps.YLim();
  //Range1D_t x = kps.XLim();

  const double logxmin = TMath::Log10(1E-5);
  const double logxmax = TMath::Log10(x.max);
  const double logymin = TMath::Log10(y.min);
  const double logymax = TMath::Log10(y.max);
/*
  double dy=0;
  double log10Ev = TMath::Log10(Ev);
  double yc = TMath::Power(10,-0.5813-0.8492*log10Ev);
  const double logymin = TMath::Log10( TMath::Max(y.min,yc-dy) );
  const double logymax = TMath::Log10( TMath::Min(y.max,yc+dy) );
*/
  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  for(int i=0; i<Nx; i++) {
   double gx = TMath::Power(10, logxmin + i * dlogx);
   for(int j=0; j<Ny; j++) {
     double gy = TMath::Power(10, logymin + j * dlogy);
     AMT = in->InitState().Tgt().Mass();  
     in->KinePtr()->Setx(gx);
     in->KinePtr()->Sety(gy);

     double xsec = fXSecModel->XSec(in, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("COHRhoKinematics", pDEBUG)  
     	 << "xsec(x= " << gx << ", y= " << gy << ") = " << xsec;
#endif
     max_xsec = TMath::Max(max_xsec, xsec);

   }//y
  }//x

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();
#endif

  return max_xsec;
}
/*void COHRhoKinematicsGenerator::dgdrot(double P[4], double costh, double sinth, double cosph, double sinph){
double P1, P2, P3;
P1=P[1];
P2=P[2];
P3=P[3];
P[1]=P1*costh*cosph-P2*sinph+P3*sinth*cosph;
P[2]=P1*costh*sinph+P2*cosph+P3*sinth*sinph;
P[3]=-P1*sinth+P3*costh;

}
*/

double COHRhoKinematicsGenerator::ComputeXYLim(const Interaction * interaction, double AMMES, Range1D_t* varx, Range1D_t* vary, double* gx, double* gy, double* gt, double* gQ2, double* gW, int* index_xy)  const
{
const InitialState & init_state = interaction->InitState();
double Ev    = init_state.ProbeE(kRfLab);
double AMT = interaction->InitState().Tgt().Mass();   //nucleus mass
double ATM=(double) init_state.Tgt().A();    //get the atm number

//std::cout<<"COHRhoKinematics/ComputeXYLim] ATM= "<<ATM<<std::endl;
double RA=R0*pow(ATM, 1.0/3.0);
double BSLOPE=(1.0/3.0)*pow((RA*FMTOGEV), 2.0);


double AMT2=AMT*AMT;

//double ATM2=ATM*ATM;           
double gy_temp=0;
double gx_temp=0;
double gt_temp=0;

Range1D_t varx_temp;
Range1D_t vary_temp;
double ymax, ymin, xmax;
const double xmin = kASmallNum;
double Elep, Amlep, Amlep2, Pklep;
double S0, Q2, Q, UU, W2,W;
RandomGen * rnd = RandomGen::Instance();

//bool nonacceptxy=true;
unsigned int Nthtry=0;
int ibadev;
do{

ibadev=0;

Nthtry=Nthtry+1;
//generate and temp gy using a random number
gy_temp = rnd->RndKine().Rndm();
gx_temp = rnd->RndKine().Rndm();

double AMMES2=AMMES*AMMES;
S0=(AMT+AMMES)*(AMT+AMMES);
xmax=(AMT2+2.0*AMT*Ev*gy_temp-S0)/(2.0*AMT*Ev*gy_temp);
if(gx_temp>xmax) {ibadev=1;}


Elep=Ev*(1.0-gy_temp);

Amlep=kMuonMass;
Amlep2=Amlep*Amlep;

if(Elep<Amlep){ibadev=2;}

Pklep=TMath::Sqrt(Elep*Elep-Amlep2);

double AA=1.0-Pklep/Elep;
double BB=1.0+Pklep/Elep;
ymax=(Ev*BB-Amlep2/(2.0*Ev))/(Ev*BB+AMT*gx_temp);
ymin=(Ev*AA-Amlep2/(2.0*Ev))/(Ev*AA+AMT*gx_temp);

if(gy_temp>ymax){ibadev=3;}
if(gy_temp<ymin){ibadev=4;}

UU=Ev*gy_temp;
Q2=2.0*AMT*UU*gx_temp;
Q=TMath::Sqrt(UU*UU+Q2);

W2=2.0*AMT*UU-Q2+AMT2;
W=TMath::Sqrt(W2);
if(W2<=S0){ibadev=5;}
if(Pklep<=0 || Q<=0){
ibadev=6;
}

double EN2C=(W2+AMT2-AMMES2)/(2.0*W);
double PN2C=sqrt(EN2C*EN2C-AMT2);

double EN1C=(W2+AMT2+Q2)/(2.0*W);
double PN1C=sqrt(EN1C*EN1C-AMT2);

double TMIN=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C-PN1C)*(PN2C-PN1C);
double TMAX=-(EN2C-EN1C)*(EN2C-EN1C)+(PN2C+PN1C)*(PN2C+PN1C);

//double TMIN2=((AMMES2+Q2)/(2.0*UU))*((AMMES2+Q2)/(2.0*UU));
gt_temp = - (1.0/BSLOPE)*log((double) rand() / double (RAND_MAX));
if(gt_temp<TMIN || gt_temp>TMAX){
//ibadev=10;
}

//make sure gy_temp is with in the range between ymax and ymin
//if not, then regenerate gy_temp
}while(ibadev !=0 && Nthtry<kRjMaxIterations);


*gy=gy_temp;
*gx=gx_temp;
varx_temp.min=xmin;
vary_temp.min=ymin;
varx_temp.max=xmax;
vary_temp.max=ymax;
*varx=varx_temp;
*vary=vary_temp;
*index_xy=Nthtry;
*gt=gt_temp;
*gQ2=Q2;
*gW=W;
//finally, get the gx, gy, ymin, ymax, xmin, xmax

return 0;
}

//___________________________________________________________________________
double COHRhoKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}

void COHRhoKinematicsGenerator::throwOnTooManyIterations(unsigned int iters,
                                                      GHepRecord* evrec) const
{
  LOG("COHKinematics", pWARN)
    << "*** Could not select valid kinematics after "
    << iters << " iterations";
  evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("Couldn't select kinematics");
  exception.SwitchOnFastForward();
  throw exception;
}



//___________________________________________________________________________
void COHRhoKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHRhoKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHRhoKinematicsGenerator::LoadConfig(void)
{
  //GetParam( "COH-Ro", fRo );
  //-- max xsec safety factor (for rejection method) and min cached energy
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6 ) ;
  // min cached energy
  GetParamDef( "Cache-MinEnergy", fEMin,  -1.0 ) ;

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);


  //-- Envelope employed when importance sampling is used 
  //   (initialize with dummy range)
  if(fEnvelope) delete fEnvelope;
  fEnvelope = new TF2("envelope",
    	  kinematics::COHImportanceSamplingEnvelope,0.,1,0.,1,2);
  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fEnvelope);
}
//____________________________________________________________________________

