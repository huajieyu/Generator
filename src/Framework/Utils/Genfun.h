//____________________________________________________________________________
/*!

\namespace genie::utils::bwfunc

\brief     Breit Wigner functions

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

\created   November 22, 2004

\cpright   Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENFUN_UTILS_H_
#define _GENFUN_UTILS_H_

namespace genie  {
namespace utils  {
namespace genfun {

  //-- A realistic Breit-Wigner distribution with L-dependent width.
  double xgencal(double XMI, double XMA,  double *XGEN);

} // genfun namespace
} // utils  namespace
} // genie  namespace

#endif   // _GENFUN_UTILS_H_
