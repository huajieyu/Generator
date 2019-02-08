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
   Moved into the new Coherent package from its previous location (EVGModules 
   package)
 @ Mar 03, 2009 - CA
   Renamed COHPiHadronicSystemGenerator -> COHRhoHadronicSystemGenerator in
   anticipation of reusing the code for simulating coherent production of
   vector mesons.
 @ Apr 02, 2009 - CA,HG,PK
   Bug fix: Reverse the order of the pion momentum rotations: Randomize the
   transverse component direction in the x'y' plane before aligning z' with 
   the direction of the momentum transfer q in the LAB.
*/
//____________________________________________________________________________

#include <cstdlib>
#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TVector3.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/Coherent/EventGen/COHRhoHadronicSystemGenerator.h"
#include "Physics/Decay/DecayModelI.h"
#include "Physics/Coherent/EventGen/COHRhoDecayer.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
COHRhoHadronicSystemGenerator::COHRhoHadronicSystemGenerator() :
HadronicSystemGenerator("genie::COHRhoHadronicSystemGenerator")
{

}
//___________________________________________________________________________
COHRhoHadronicSystemGenerator::COHRhoHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::COHRhoHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
COHRhoHadronicSystemGenerator::~COHRhoHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void COHRhoHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system (rho + nucleus) in 
// COH interactions
//
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();

  if (fXSecModel->Id().Name() == "genie::VectDominCOHRhoXSec") {
    CalculateHadronicSystem_VectDomin(evrec);
  }    
  // Add the coherent rho decay products
  //int pdgc = getRhoPDGCode(evrec);
  //this->AddRhoDecayProducts(evrec, pdgc);

  //LOG("COHRhoHadronicVtx", pNOTICE)
  //   << "Decay coherent rho in  the initial hadronic decay products";
  //this->PreHadronTransportDecays(evrec);

/*
*/
} 
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void COHRhoHadronicSystemGenerator::CalculateHadronicSystem_VectDomin(GHepRecord * evrec) const
{
  // This method generates the final state hadronic system (pion + nucleus) in 
  // COH interactions
  //   
  RandomGen * rnd = RandomGen::Instance();

  Interaction * interaction = evrec->Summary();
  const XclsTag & xcls_tag  = interaction->ExclTag();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  //const TLorentzVector & vtx   = *(nu->X4());
  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());
  
  //-- Determine the pdg code of the final state pion & nucleus
  int nucl_pdgc = Ni->Pdg(); // same as the initial nucleus
  //int pion_pdgc = getPionPDGCodeFromXclTag(xcls_tag);
      
  //basic kinematic inputs
  double E    = nu->E();
  //double M    = kNucleonMass;
  //double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  //double mpi2 = TMath::Power(mpi,2);
  double xo   = interaction->Kine().x(true);
  double yo   = interaction->Kine().y(true);
  double to   = interaction->Kine().t(true);
  double Q2   = interaction->Kine().Q2(true);
  double AMT = interaction->InitState().Tgt().Mass();   //nucleus mass
  //double AMT2 = AMT*AMT;
  double AMMES = interaction->ExclTag().RhoMass();
  double AMMES2 = AMMES*AMMES;
  double v=E*yo;
  double Q=TMath::Sqrt(v*v+Q2);


  SLOG("COHHadronicVtx", pINFO)
    << "Ev = "<< E << ", xo = " << xo
    << ", yo = " << yo << ", to = " << to;

  // 4-momentum transfer q=p(neutrino)-p(f/s lepton)
  TLorentzVector q=p4nu-p4fsl;
  SLOG("COHHadronicVtx", pINFO)
    << "\n 4-p transfer q @ LAB: " << utils::print::P4AsString(&q);


  double PNf[4], PMES[4];
  PNf[0]=AMT+to/(2.0*AMT);
  //double PKNf=TMath::Sqrt(PNf[0]*PNf[0]-AMT2);

  double PKMES, PLMES, PTMES;
  PMES[0]=E*yo+AMT-PNf[0];
  PKMES=TMath::Sqrt(PMES[0]*PMES[0]-AMMES2);
  double sinx1, cosx1,sinph, cosph, sinla, cosla;
  cosx1=(Q2-to-AMMES2+2*PMES[0]*v)/(2.0*Q*PKMES);
  if(TMath::Abs(cosx1)>1.0) {cosx1=1.0; sinx1=0.0;}
  sinx1=TMath::Sqrt(1.0-cosx1*cosx1);

  PLMES=PKMES*cosx1;
  PTMES=PKMES*sinx1;
  double phi=2.0*kPi*rnd->RndHadro().Rndm(); 
  sinph=sin(phi); cosph=cos(phi);

  PMES[1]=PTMES*sinph;
  PMES[2]=PTMES*cosph;
  PMES[3]=PLMES;

//  double Q=TMath::Sqrt(v*v+Q2);

  double PLNUC=Q-PLMES;
  double PTNUC=-PTMES;

  PNf[1]=PTNUC*sinph;
  PNf[2]=PTNUC*cosph;
  PNf[3]=PLNUC;

  double PTLEP=TMath::Sqrt(p4fsl.Px()*p4fsl.Px()+p4fsl.Py()*p4fsl.Py());
  double sinal=-PTLEP/Q;
  double cosal=0;
  if(TMath::Abs(sinal)>1.0) {sinal=1.0; cosal=0.0;}
  cosal=TMath::Sqrt(1.0-sinal*sinal);

  //double LAM=2.0*kPi*((double) rand() / double (RAND_MAX));
  //cosla=cos(LAM);
  //sinla=sin(LAM);
  //get sinla and cosla from the momentum of the leading lepton
  cosla=p4fsl.Px()/PTLEP; 
  sinla=p4fsl.Py()/PTLEP;


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

  TLorentzVector PNff, PMESf;
  //PNff.SetPx(PNf[1]);
  //PNff.SetPy(PNf[2]);
  //PNff.SetPz(PNf[3]);
  //PNff.SetE(PNf[0]);

  //PMESf.SetPx(PMES[1]);
  //PMESf.SetPy(PMES[2]);
  //PMESf.SetPz(PMES[3]);
  //PMESf.SetE(PMES[0]);
  PNff.SetPxPyPzE(PNf[1],PNf[2], PNf[3], PNf[0]);
  PMESf.SetPxPyPzE(PMES[1],PMES[2],PMES[3],PMES[0]);



  //--Save the particles at the GHEP record;
  int mom = evrec->TargetNucleusPosition();
  int rho_pdgc = interaction->ExclTag().RhoPDG();        
  evrec->AddParticle(
         rho_pdgc,kIStStableFinalState, mom,-1,-1,-1,
         PMESf.Px(), PMESf.Py(), PMESf.Pz(), PMESf.E(), 0, 0, 0, 0); 
  
  
  //int nucl_pdgc = interaction->InitState().TgtPdg(); // same as the initial nucleus 
  evrec->AddParticle(
        nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1,
        PNff.Px(), PNff.Py(), PNff.Pz(), PNff.E(), 0, 0, 0, 0); 

   //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
   




}
///////////////////////////////////////////////////////////
int COHRhoHadronicSystemGenerator::getPionPDGCodeFromXclTag(const XclsTag& xcls_tag) const
{
  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     == 1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  == 1) pion_pdgc = kPdgPiP;
  else if (xcls_tag.NPiMinus() == 1) pion_pdgc = kPdgPiM;
  else {
    LOG("COHHadronicVtx", pFATAL)
      << "No final state pion information in XclsTag!";
    exit(1);
  }
  return pion_pdgc;
}
//___________________________________________________________________________
int COHRhoHadronicSystemGenerator::getRhoPDGCode(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  
  // get rho id
  const XclsTag & xls = interaction->ExclTag(); 
  int pdgc=xls.RhoPDG();

  return pdgc;
}

//____________________________________________________________________

//_________________________________________________________________________
void COHRhoHadronicSystemGenerator::AddRhoDecayProducts(
                                         GHepRecord * evrec, int pdgc) const
{
  //Decay the rho state, take the decay products and 
  //boost them into LAB frame;
  // find the rho position
  int irpos = evrec->ParticlePosition(pdgc, kIStDecayedState, 0);
  assert(irpos>0);

  // access the GHEP entry
  GHepParticle * cohrho = evrec->Particle(irpos);
  assert(cohrho);

  // resonance location
  const TLorentzVector & x4 = *(cohrho->X4());

  const TLorentzVector & p4 = *(cohrho->P4());
  
  //prepare the decayer inputs
  DecayerInputs_t dinp;
  dinp.PdgCode=pdgc;

  dinp.P4 =cohrho->P4();

  //do the decay
  
  TClonesArray * decay_products = fRhoDecayer->Decay(dinp);
  
    if(!decay_products) {
     LOG("RESHadronicVtx", pWARN) << "Got an empty decay product list!";
     LOG("RESHadronicVtx", pWARN)
                      << "Quitting the current event generation thread";

     evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);

     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Not enough phase space for hadronizer");
     exception.SwitchOnFastForward();
     throw exception;

     return;
    }
   
    // decide the istatus of decay products
    GHepParticle * nuc = evrec->TargetNucleus(); 
    GHepStatus_t dpist = (nuc) ? kIStHadronInTheNucleus : kIStStableFinalState;
    
    if(decay_products) 
    {
     cohrho->SetStatus(kIStDecayedState);
     // loop over the daughter and add them to the event record
     TMCParticle * dpmc =0;
     TObjArrayIter decay_iter(decay_products);
     while ( (dpmc = (TMCParticle *) decay_iter.Next()) ) {
     int dppdg = dpmc->GetKF();
     double px = dpmc->GetPx();
     double py = dpmc->GetPy();
     double pz = dpmc->GetPz();
     double E  = dpmc->GetEnergy();
     TLorentzVector p4(px,py,pz,E); 
     
     // only add the decay products, the mother particle already exists
     if(dpmc->GetKS()==1) {
         evrec->AddParticle(dppdg,dpist,irpos,-1,-1,-1, p4, x4);
     }
     }
    //done release the original list
    decay_products->Delete();
    delete decay_products;
    }
    
   //--------------------------------------------------------------------
    double rho_px, rho_py, rho_pz, rho_E, AMMES;
    double Epi, Ppi1[4], Ppi0[4],Pkpi;
    double DELTA2, DELTA1;
    rho_px=p4.Px();
    rho_py=p4.Py();
    rho_pz=p4.Pz();
    rho_E=p4.E();
    AMMES=TMath::Sqrt(rho_E*rho_E-rho_px*rho_px-rho_py*rho_py-rho_pz*rho_pz);
    Epi=AMMES/2.0;
    Pkpi=TMath::Sqrt(Epi*Epi-kPionMass*kPionMass);
    DELTA2=2.0*kPi*((double) rand() / double (RAND_MAX));
    DELTA1=acos(((double) rand() / double (RAND_MAX))*2.0-1.0);
       
    Ppi1[1]=Pkpi*sin(DELTA1)*cos(DELTA2);
    Ppi1[2]=Pkpi*sin(DELTA1)*sin(DELTA2);
    Ppi1[3]=Pkpi*cos(DELTA1);
    Ppi1[0]=Epi;

    Ppi0[1]=-Ppi1[1];
    Ppi0[2]=-Ppi1[2];
    Ppi0[3]=-Ppi1[3];
    Ppi0[0]=Epi; 
    //lorentz transform to lab frame
    TLorentzVector v_pip, v_pi0, CMS;

    v_pip.SetPxPyPzE(Ppi1[1], Ppi1[2], Ppi1[3], Epi);
    v_pi0.SetPxPyPzE(Ppi0[1], Ppi0[2], Ppi0[3], Epi);
    CMS=p4;

    v_pip.Boost(CMS.BoostVector());
    v_pi0.Boost(CMS.BoostVector());

    //reverset the direction of the pions
    Ppi1[1]=v_pip.Px();
    Ppi1[2]=v_pip.Py();
    Ppi1[3]=v_pip.Pz();

    Ppi0[1]=v_pi0.Px();
    Ppi0[2]=v_pi0.Py();
    Ppi0[3]=v_pi0.Pz();

    v_pip.SetPxPyPzE(Ppi1[1], Ppi1[2], Ppi1[3], Epi);
    v_pi0.SetPxPyPzE(Ppi0[1], Ppi0[2], Ppi0[3], Epi);

     
    int dppdg[2];
    dppdg[0]=111;
    dppdg[1]=211;
   // evrec->AddParticle(dppdg[0],dpist,irpos,-1,-1,-1, v_pi0, x4);
   // evrec->AddParticle(dppdg[1],dpist,irpos,-1,-1,-1, v_pip, x4);
    
   //-------------------------------------------------------------------e

}
//_________________________________________________________________________

void COHRhoHadronicSystemGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________
void COHRhoHadronicSystemGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void COHRhoHadronicSystemGenerator::LoadConfig(void)
{
fRhoDecayer = 0;
fRhoDecayer =
         dynamic_cast<const DecayModelI *> (this->SubAlg("Decayer"));
assert(fRhoDecayer);

}


