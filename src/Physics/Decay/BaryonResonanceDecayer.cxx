//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 27, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 02, 2009 - CA
   Add dummy `UnInhibitDecay(int,TDecayChannel*) const' and `InhibitDecay(int,
   TDecayChannel*) const' methods to conform to the DecayModelI interface.
   To implement soon.
*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TMath.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/Decay/BaryonResonanceDecayer.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer() :
DecayModelI("genie::BaryonResonanceDecayer")
{

}
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer(string config) :
DecayModelI("genie::BaryonResonanceDecayer", config)
{

}
//____________________________________________________________________________
BaryonResonanceDecayer::~BaryonResonanceDecayer()
{

}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsHandled(int code) const
{
// handles only requests to decay baryon resonances

  if( utils::res::IsBaryonResonance(code) ) return true;

  LOG("Decay", pINFO) 
      << "This algorithm can not decay particles with PDG code = " << code;

  return false;
}
//____________________________________________________________________________

//-----------------------------------------customizing-----------------------
double BaryonResonanceDecayer::BRDeltaNGammaNPi(int id_mother, int ichannel, double W) const
{
  //-- auxiliary parameters
  int DeltaFlag = 0;
  if (id_mother == 2114 || id_mother==-2114) {
    DeltaFlag = 1; // Delta0 or Delta0_bar
  }  
  else if (id_mother == 2214 || id_mother==-2214) {
    DeltaFlag = 2; // Delta+ or anti_Delta+
  }  
  else  {
    return 0;
  }

  if( W <= genie::constants::kNucleonMass + genie::constants::kPi0Mass ) {
    if (ichannel == 0) {return 0;} // ichannel =0,1,2 has to match
    // the channel order in genie_pdg_table.dat
    if (ichannel == 1) {return 0;}
    if (ichannel == 2) {return 1;}
  } else {

  // double rDelta= 0.81*FMTOGEV;

    double W_2   = TMath::Power(W,    2);

    double pPiW   = TMath::Sqrt( (W_2 - fMAux1) * (W_2 - fMAux2) )/(2*W);
    double pPim   = TMath::Sqrt( ( fDeltaMass2- fMAux1 )*(fDeltaMass2- fMAux2) )/( 2 * fDeltaMass );

    double EgammaW = (W_2-fNucleonMass2)/(2*W);
    double Egammam = (fDeltaMass2-fNucleonMass2)/(2*fDeltaMass);

    double TPiW=TMath::Power(pPiW, 3);
    double TPim=TMath::Power(pPim, 3);

    double fgammaW= 1./(TMath::Power(1+EgammaW*EgammaW/fFFScaling, 2));
    double fgammam= 1./(TMath::Power(1+Egammam*Egammam/fFFScaling, 2));

    double Rinverse = fWidthPi_0 * TMath::Power(Egammam, 3)*TMath::Power(fgammam, 2)*TPiW
                         /(fWidthGamma_0*TMath::Power(EgammaW, 3)*TMath::Power(fgammaW, 2)*TPim);

    double BRPi = Rinverse/(1+Rinverse);
    double BRgamma = 1/(1+Rinverse);
  
    if (DeltaFlag==1) {
      if (ichannel == 0) {return BRPi * fBRPi02 ;}
      if (ichannel == 1) {return BRPi * fBRPi01 ;}
      if (ichannel == 2) {return BRgamma;}
    }
    if (DeltaFlag==2) {
      if (ichannel == 0) {return BRPi * fBRPi01 ;}
      if (ichannel == 1) {return BRPi * fBRPi02 ;}
      if (ichannel == 2) {return BRgamma;}
    }
  }

  return 0;

}
//-----------------------------------------------
TClonesArray* BaryonResonanceDecayer::Decay(const DecayerInputs_t & inp) const
{  
  if ( ! this->IsHandled(inp.PdgCode) ) return 0;

  //-- Find the particle in the PDG library & quit if it does not exist
  TParticlePDG * mother = PDGLibrary::Instance()->Find(inp.PdgCode);

  if(!mother) {
     LOG("Decay", pERROR)
          << "\n *** The particle with PDG-Code = " << inp.PdgCode
                                         << " was not found in PDGLibrary";
     return 0;                               
  }  
  LOG("Decay", pINFO)
       << "Decaying a " << mother->GetName()
                        << " with P4 = " << utils::print::P4AsString(inp.P4);
  
  //-- Reset previous weight
  fWeight = 1.;

  //-- Get the resonance mass W (generally different from the mass associated
  //   with the input pdg_code, since the it is produced off the mass shell)
  double W = inp.P4->M();
  LOG("Decay", pINFO) << "Available mass W = " << W;
  
  //-- Get all decay channels
  TObjArray * decay_list = mother->DecayList();
  unsigned int nch = decay_list->GetEntries();
  LOG("Decay", pINFO)
               << mother->GetName() << " has: " << nch << " decay channels";

  //-- Loop over the decay channels (dc) and write down the branching
  //   ratios to be used for selecting a decay channel.
  //   Since a baryon resonance can be created at W < Mres, explicitly
  //   check and inhibit decay channels for which W > final-state-mass
  
  double BR[nch], tot_BR = 0;    
 
  //------------------ cusomizing ------------------------------
  if(inp.PdgCode== 2114 || inp.PdgCode==-2114 || inp.PdgCode==2214||inp.PdgCode==-2214){
	  for(unsigned int ich = 0; ich < nch; ich++) {

     TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
     double fsmass = this->FinalStateMass(ch);

     if(fsmass < W) {
       SLOG("Decay", pDEBUG)
               << "Using channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
       tot_BR += BaryonResonanceDecayer::BRDeltaNGammaNPi(inp.PdgCode, ich, W);
     } else {       
       SLOG("Decay", pINFO)
               << "Suppresing channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
     }
     BR[ich] = tot_BR;
  }
  }  else {
  for(unsigned int ich = 0; ich < nch; ich++) {

     TDecayChannel * ch = (TDecayChannel *) decay_list->At(ich);
     double fsmass = this->FinalStateMass(ch);

     if(fsmass < W) {
       SLOG("Decay", pDEBUG)
               << "Using channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
       tot_BR += ch->BranchingRatio();
     } else {       
       SLOG("Decay", pINFO)
               << "Suppresing channel: " << ich 
                        << " with final state mass = " << fsmass << " GeV";         
     }
     BR[ich] = tot_BR;
  }
  }
//--------------customizing ends --------------------------------------------------
/*  
*/
  if(tot_BR==0) {
    SLOG("Decay", pWARN) 
      << "None of the " << nch << " decay chans is available @ W = " << W;
    return 0;    
  }

  //-- Select a resonance based on the branching ratios
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = tot_BR * rnd->RndDec().Rndm();
  do { 
    sel_ich = ich;
    
  } while (x > BR[ich++]);

  TDecayChannel * ch = (TDecayChannel *) decay_list->At(sel_ich);

  LOG("Decay", pINFO) 
    << "Selected " << ch->NDaughters() << "-particle decay chan (" 
    << sel_ich << ") has BR = " << ch->BranchingRatio();

  //-- Decay the exclusive state and return the particle list
  TLorentzVector p4(*inp.P4);

  double Q2=inp.Qsqr;
  TLorentzVector p4v(*inp.P4v);
  TLorentzVector p4l(*inp.P4l);

  return ( this->DecayExclusive(inp.PdgCode, p4, ch, Q2, p4v, p4l) );
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * BaryonResonanceDecayer::DecayExclusive(
     int pdg_code, TLorentzVector & p, TDecayChannel * ch,double Q2, TLorentzVector p4v, TLorentzVector p4l) const
{
  //-- Get the final state mass spectrum and the particle codes
  unsigned int nd = ch->NDaughters();

  int    pdgc[nd];
  double mass[nd];

  bool DeltaNpidecay=false;      // flag of expected channel Delta->pion+nucleon
  int  npi=0,  nnucl=0;
  unsigned int pi_id = 0;

  //the propability of decaying into three body is very small, simulated in GENIE? 
  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     pdgc[iparticle] = daughter_code;
     mass[iparticle] = daughter->Mass();

     if ( daughter_code == 211 || daughter_code == 111 || daughter_code == -211 ) {
       ++npi ;
       pi_id = iparticle ;
     }
     if ( daughter_code == 2112 || daughter_code == 2212 ) {
       ++nnucl ;
     }

     SLOG("Decay", pINFO)
       << "+ daughter[" << iparticle << "]: "
        << daughter->GetName() << " (pdg-code = "
          << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
  }  // loop over the duaghter particles

  if( nd==2 && (pdg_code==2224 || pdg_code==2214 || pdg_code==2114) ) {
    if( npi==1 && nnucl==1 ){
      DeltaNpidecay=true;
    }
  }

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(p, nd, mass);
  assert(is_permitted);

  // ANL OR BNL ANGULAR DISTRIBUTION(source comes from NuWRO)-----------------------------------------------------
  double r33,r31,r3m1;
  double ct,cphi,weight;
  double maxwei=1.;
  
  switch(fANLORBNL)
  {
     //No correlations
     case 0:
     {
          r33=0.5;
          r31=0.;
          r3m1=0.;
          maxwei=1.;
          break;
     }

     //ANL
     case 1:
     {
       if(Q2<0.1)
       {
          r33=0.523;
          r31=-0.322;
          r3m1=-0.138;
          maxwei=1.+(r33-0.5)-2.*sqrt(3.)*r31-sqrt(3.)*r3m1;

          break;
       } 
       else if(Q2<0.3)
       {
          r33=0.649;
          r31=-0.128;
          r3m1=0.034;
          maxwei=1.+(r33-0.5)-2.*sqrt(3.)*r31+sqrt(3.)*r3m1;

          break;
       } 
       else if(Q2<0.5)
       {
          r33=0.674;
          r31=-0.017;
          r3m1=-0.203;
          maxwei=1.+(r33-0.5)-2.*sqrt(3.)*r31-sqrt(3.)*r3m1;

          break;
       } 
       else
       {
          r33=0.748;
          r31=0.041;
          r3m1=-0.162;
          maxwei=1.+(r33-0.5)+2.*sqrt(3.)*r31-sqrt(3.)*r3m1;

          break;
       } 
     }   

     //BNL
     case 2:
     {
       if(Q2<0.2)
       {
          r33=0.661;
          r31=-0.213;
          r3m1=-0.133;
          maxwei=1.+(r33-0.5)-2.*sqrt(3.)*r31-sqrt(3.)*r3m1;

          break;
       } 
       else if(Q2<0.4)
       {
          r33=0.673;
          r31=0.025;
          r3m1=-0.075;
          maxwei=1.+(r33-0.5)+2.*sqrt(3.)*r31-sqrt(3.)*r3m1;

          break;
       } 
       else if(Q2<0.6)
       {
          r33=0.750;
          r31=0.036;
          r3m1=0.046;
          maxwei=1.+(r33-0.5)+2.*sqrt(3.)*r31+sqrt(3.)*r3m1;

          break;
       } 
       else
       {
          r33=0.800;
          r31=0.075;
          r3m1=0.004;
          maxwei=1.+(r33-0.5)+2.*sqrt(3.)*r31+sqrt(3.)*r3m1;

          break;
       } 
     }   
     //Theta only
     case 3:
     {
         r33=0.75;
         r31=0.;
         r3m1=0.;
         maxwei=1.25;

         break;
     }
     //No correlations
     default:
     {
          r33=0.5;
          r31=0.;
          r3m1=0.;
          maxwei=1.;

          break;
     }

  }

  TLorentzVector q=p4v-p4l;


  //-- Create the event record
   TClonesArray * particle_list = new TClonesArray("TMCParticle", 1+nd) ;

   //-- Add the mother particle to the event record (KS=11 as in PYTHIA)
   TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);

   new ( (*particle_list)[0] )
        TMCParticle( 11,pdg_code,
                     0,0,0,
                     p.Px(), p.Py(), p.Pz(), p.Energy(), mother -> Mass(),
                     0,0,0,0,0 ) ;

  double wmax = -1, w ;
  for( int i=0; i<50; i++ ) {
     w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("Decay", pINFO)
     << "Max phase space gen. weight for current decay: " << wmax;

  if ( fGenerateWeighted )  {
     // *** generating weighted decays ***
     w = fPhaseSpaceGenerator.Generate();
     fWeight *= TMath::Max(w/wmax, 1.);
  }
  else  {
     // *** generating un-weighted decays ***
     RandomGen * rnd = RandomGen::Instance();
     wmax *= 2;
     bool accept_decay=false;
     unsigned int itry=0;

     while(!accept_decay) {
       itry++;
       assert(itry<kMaxUnweightDecayIterations);

       w  = fPhaseSpaceGenerator.Generate();
       double gw = wmax * rnd->RndDec().Rndm();

       if( w > wmax ) {
          LOG("Decay", pWARN) 
             << "Current decay weight = " << w << " > wmax = " << wmax;
       }
       LOG("Decay", pINFO) 
          << "Current decay weight = " << w << " / R = " << gw;

       accept_decay = (gw<=w);

       if ( accept_decay && DeltaNpidecay ) {
         // in this case we need additional checks based on the kinamatics of the pion and of the nucleon

         TLorentzVector * temp ;

         // retrieve the pion
         temp = fPhaseSpaceGenerator.GetDecay(pi_id);
         TLorentzVector pion( temp -> Px(), temp -> Py(), temp -> Pz(), temp -> Energy() ) ;
 
         pion.Boost(-p.BoostVector() );
          p4v.Boost(-p.BoostVector() );
          p4l.Boost(-p.BoostVector() );
            q.Boost(-p.BoostVector() );

         TVector3 p3v(p4v.Px(),p4v.Py(),p4v.Pz());
         TVector3 p3l(p4l.Px(),p4l.Py(),p4l.Pz()); 
         TVector3 mtr(q.Px(),q.Py(),q.Pz()); 

         TVector3 zaxis=1./sqrt(mtr.X()*mtr.X()+mtr.Y()*mtr.Y()+mtr.Z()*mtr.Z())*mtr;
 
         //TVector3 zaxis=mtr/sqrt(mtr.X()*mtr.X()+mtr.Y()*mtr.Y()+mtr.Z()*mtr.Z());
         TVector3 yaxis=p3v.Cross(p3l);
                  yaxis=1./sqrt(yaxis.X()*yaxis.X()+yaxis.Y()*yaxis.Y()+yaxis.Z()*yaxis.Z())*yaxis;
         //           yaxis=yaxis/sqrt(yaxis.X()*yaxis.X()+yaxis.Y()*yaxis.Y()+yaxis.Z()*yaxis.Z());
         TVector3 xaxis=yaxis.Cross(zaxis);

         //-----------------------------------------------------------------
         TVector3 pion3m(pion.Px(),pion.Py(),pion.Pz());
         TVector3 pion_dir=1./sqrt(pion.Px()*pion.Px()+pion.Py()*pion.Py()
                           +pion.Pz()*pion.Pz())*pion3m;

         ct=pion_dir*zaxis;
         cphi=pion_dir*xaxis;

         weight=(1.- (r33-0.5)*(3.*ct*ct-1.)-sqrt(3.0)*( 2.*r31*sqrt(1.-ct*ct)*ct*cphi + r3m1*(1.-ct*ct)*(2.*cphi*cphi-1.) ));
         double aidrnd=maxwei*gRandom->Rndm();

         if ( weight < aidrnd )
           accept_decay = false ;

       }  // It is a Delta -> N Pi decay

     }  // loop over possible decays

  } // not weighted generation


  // we fill the event record as it's coming from the phase space decayer

  for(unsigned int id = 0; id < nd; id++) {

    TLorentzVector * p4 = fPhaseSpaceGenerator.GetDecay(id);

    LOG("Decay", pDEBUG)
    << "Adding final state particle PDGC = " << pdgc[id]
                                                     << " with mass = " << mass[id] << " GeV";
    new ( (*particle_list)[1+id] )
              TMCParticle(1,pdgc[id],
                          0,0,0,
                          p4 -> Px(), p4 -> Py(), p4 -> Pz(), p4 -> Energy(), mass[id],
                          0,0,0,0,0 ) ;

  } //end daughter particle loop

  //-- Set owner and return
  particle_list->SetOwner(true);
  return particle_list;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

}
//____________________________________________________________________________
void BaryonResonanceDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;

  if(!dc) return;

}
//____________________________________________________________________________
double BaryonResonanceDecayer::FinalStateMass(TDecayChannel * ch) const
{
// Computes the total mass of the final state system

  double mass = 0;
  unsigned int nd = ch->NDaughters();

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

    int daughter_code = ch->DaughterPdgCode(iparticle);
    TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
    assert(daughter);

    double md = daughter->Mass();

    // hack to switch off channels giving rare  occurences of |1114| that has
    // no decay channels in the pdg table (08/2007)
    if(TMath::Abs(daughter_code) == 1114) {
      LOG("Decay", pNOTICE)
                      << "Disabling decay channel containing resonance 1114";;
      md = 999999999;
    }
    mass += md;
  }  
  return mass;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BaryonResonanceDecayer::LoadConfig(void)
{
// Read configuration options or set defaults

  //-- Generated weighted or un-weighted hadronic systems

  // note that this variable is not present in any of the xml configuration files
  fGenerateWeighted = false ;
  //GetParam( "generate-weighted", fGenerateWeighted, false );  decomment this line if the variable needs to be taken from configurations

  this->GetParam( "Prob32", fProb32 ) ;
  fProb12 = 1. - fProb32 ;

  double delta_width ;
  this -> GetParam( "DeltaWidth", delta_width ) ;

  this -> GetParam( "DeltaMass", fDeltaMass ) ;

  double TotGammaBR ;
  this->GetParam( "TotGammaBR", TotGammaBR ) ;

  double TotNPiBR = 1. - TotGammaBR ;

  fWidthPi_0 =    delta_width * TotNPiBR ;
  fWidthGamma_0 = delta_width * TotGammaBR ;

  this -> GetParam( "OnePiBR", fBRPi01 ) ;
  fBRPi02 = 1. - fBRPi01 ;

  fDeltaMass2   = TMath::Power( fDeltaMass, 2) ;
  fNucleonMass2 = TMath::Power( genie::constants::kNucleonMass, 2) ;
  fMAux1 = TMath::Power( genie::constants::kNucleonMass + genie::constants::kPi0Mass , 2) ;
  fMAux2 = TMath::Power( genie::constants::kNucleonMass - genie::constants::kPi0Mass , 2) ;

  this -> GetParam( "FFScaling", fFFScaling ) ;

  this -> GetParam("ANLORBNL", fANLORBNL);


}
//____________________________________________________________________________
