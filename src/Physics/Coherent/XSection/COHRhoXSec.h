//____________________________________________________________________________
/*!

\class    genie::COHRhoXSec

\brief    Computes the cross section for COH neutrino-nucleus pi production.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHRho_XSEC_H_
#define _COHRho_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

  class COHRhoXSec : public XSecIntegratorI {

    public:
      COHRhoXSec();
      COHRhoXSec(string config);
      virtual ~COHRhoXSec();

      // XSecIntegratorI interface implementation
      double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

      // Overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      void LoadConfig (void);
      double ComputeXYLim (const Interaction * in, double AMMES, Range1D_t* varx, Range1D_t* vary, double* gx, double* gy) const;
      double fQ2Min;  ///< lower bound of integration for Q^2 in Berger-Sehgal Model
      double fQ2Max;  ///< upper bound of integration for Q^2 in Berger-Sehgal Model
      double fTMax;   ///< upper bound for t = (q - p_pi)^2
  };

}       // genie namespace
#endif  // _COHRho_XSEC_H_
