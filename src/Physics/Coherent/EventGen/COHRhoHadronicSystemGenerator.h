//____________________________________________________________________________
/*!

\class    genie::COHRhoHadronicSystemGenerator

\brief    Generates the f/s hadronic system in v COH pi production interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHRHO_HADRONIC_SYSTEM_GENERATOR_H_
#define _COHRHO_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

  class XclsTag;
  class DecayModelI;
  class COHRhoHadronicSystemGenerator : public HadronicSystemGenerator {

  public :
    COHRhoHadronicSystemGenerator();
    COHRhoHadronicSystemGenerator(string config);
    ~COHRhoHadronicSystemGenerator();

    // implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const;
    void CalculateHadronicSystem_VectDomin(GHepRecord * event_rec) const;
    void Configure(const Registry & config);
    void Configure(string config);
  private:
    void LoadConfig                (void);
    int getPionPDGCodeFromXclTag(const XclsTag& xcls_tag) const;
    int getRhoPDGCode(GHepRecord * evrec) const; 
    void AddRhoDecayProducts (GHepRecord * evrec, int pdgc) const;
    const DecayModelI * fRhoDecayer;



};

}      // genie namespace
#endif // _COHRHO_HADRONIC_SYSTEM_GENERATOR_H_

