//____________________________________________________________________________
/*!

\class    genie::COHRhoKinematicsGenerator

\brief    Generates values for the kinematic variables describing coherent 
          neutrino-nucleus pion production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHRHO_KINEMATICS_GENERATOR_H_
#define _COHRHO_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

class TF2;

namespace genie {

class COHRhoKinematicsGenerator : public KineGeneratorWithCache {

public :
  COHRhoKinematicsGenerator();
  COHRhoKinematicsGenerator(string config);
 ~COHRhoKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:
  // methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);

  // overload KineGeneratorWithCache method to compute max xsec
  double ComputeMaxXSec (const Interaction * in) const;
  //void   dgdrot(double P[4],double cosal, double sinal, double cosla, double sinla) const;
  double ComputeXYLim (const Interaction * in, double AMMES, Range1D_t* varx, Range1D_t* vary, double* gx, double* gy,double* gt,double* gQ2, double* gW, int* index_xy) const;
  // overload KineGeneratorWithCache method to get energy
  double Energy         (const Interaction * in) const;

  mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
  double fRo;              ///< nuclear scale parameter
private:
  void   throwOnTooManyIterations(unsigned int iters, GHepRecord* evrec) const; 

};

}      // genie namespace
#endif // _COHRHO_KINEMATICS_GENERATOR_H_
