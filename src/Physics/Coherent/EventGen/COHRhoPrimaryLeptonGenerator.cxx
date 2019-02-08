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

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/Coherent/EventGen/COHRhoPrimaryLeptonGenerator.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
COHRhoPrimaryLeptonGenerator::COHRhoPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::COHRhoPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
COHRhoPrimaryLeptonGenerator::COHRhoPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::COHRhoPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
COHRhoPrimaryLeptonGenerator::~COHRhoPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void COHRhoPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in COH events

//  Interaction * interaction = evrec->Summary();

 // const InitialState & init_state = interaction->InitState();
/*
  GHepParticle * mom = evrec->Probe();
  int imom = evrec->ProbePosition();


  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);
  //double x  = interaction->Kine().x(true);

  // Auxiliary params
  double Ev  = init_state.ProbeE(kRfLab);
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);
*/

/*
  double sinth=0.0; 
  double costh=0.0;
  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double pklep=TMath::Sqrt(El*El-ml*ml);
  costh=(2*Ev*El-ml2-Q2)/(2*Ev*pklep);
  if(TMath::Abs(costh)>1.0) {costh=1.0; sinth=0.0;}
  sinth=TMath::Sqrt(1.0-costh*costh);

  double plp = pklep*costh;                         // p(//)
  double plt = pklep*sinth;                         // p(-|)

  LOG("LeptonicVertex", pNOTICE)
          << "fsl: E = " << El << ", |p//| = " << plp << "[pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  //p3l.RotateUz(unit_nudir);

  
  // Lepton 4-momentum in LAB
  //TVector3 p3l(pltx, plty, plp);
  TLorentzVector p4l(p3l,El);
  TLorentzVector & vtx= *(mom->X4());
  TLorentzVector x4l(vtx);
*/
  // Figure out the Final State Lepton PDG Code
  //int pdgc = interaction->FSPrimLepton()->PdgCode();

  //evrec->AddParticle(pdgc, kIStStableFinalState, imom, -1,-1,-1, p4l, x4l);
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);//
  // Create a GHepParticle and add it to the event record
  //this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  //this->SetPolarization(evrec);

}
//___________________________________________________________________________
