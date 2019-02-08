//____________________________________________________________________________
/*!

\class    genie::COHRhoPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v COH NC interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 26, 2005

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHRHO_PRIMARY_LEPTON_GENERATOR_H_
#define _COHRHO_PRIMARY_LEPTON_GENERATOR_H_

#include "Physics/Common/PrimaryLeptonGenerator.h"

namespace genie {

class COHRhoPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  COHRhoPrimaryLeptonGenerator();
  COHRhoPrimaryLeptonGenerator(string config);
  ~COHRhoPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _COHRHO_PRIMARY_LEPTON_GENERATOR_H_
