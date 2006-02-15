//____________________________________________________________________________
/*!

\class    genie::QELXSec

\brief    Computes the Quasi Elastic (QEL) cross section.

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecAlgorithmI("genie::QELXSec")
{

}
//____________________________________________________________________________
QELXSec::QELXSec(string config) :
XSecAlgorithmI("genie::QELXSec", config)
{

}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::XSec(const Interaction * in) const
{
  if(! this -> ValidProcess    (in) ) return 0.;
  if(! this -> ValidKinematics (in) ) return 0.;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Get initial & final state information
  const InitialState & init_state = interaction->GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);

  // Estimate the integration limits & step
  Range1D_t  rQ2 = utils::kinematics::Q2Range_M(interaction);
  LOG("QELXSec", pDEBUG) << "Q2 integration range = ("
                                    << rQ2.min << ", " << rQ2.max << ")";
  double logQ2max = TMath::Log(rQ2.max);
  double logQ2min = TMath::Log(rQ2.min);
  double dlogQ2   = (logQ2max - logQ2min) / (fNBins - 1);

  // Define the integration grid & instantiate a FunctionMap
  UnifGrid grid;
  grid.AddDimension(fNBins, logQ2min, logQ2max);

  FunctionMap Q2dxsec_dQ2(grid);

  // Loop over logQ2 range, estimate/store Q2*dxsec/dQ2
  for(int iQ2 = 0; iQ2 < fNBins; iQ2++) {
     double Q2 = TMath::Exp(logQ2min + iQ2 * dlogQ2);

     //-- update the scattering parameters
     interaction->GetKinematicsPtr()->SetQ2(Q2);
     //-- compute dxsec/dQ2
     double pxsec = fDiffXSecModel->XSec(interaction);
     //-- push Q2*(dxsec/dQ2) to the FunctionMap
     Q2dxsec_dQ2.AddPoint(Q2*pxsec, iQ2);

     LOG("QELXSec", pDEBUG)
          << "point...." << iQ2+1 << "/" << fNBins << " : "
                  << "dxsec/dQ^2 (Q^2 = " << Q2 << " ) = " << pxsec;
  }
  delete interaction;

  // Do the numerical integration
  double xsec = fIntegrator->Integrate(Q2dxsec_dQ2);
  LOG("QELXSec", pDEBUG) << "XSec[QEL] (E = " << E << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool QELXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  int  nu  = init_state.GetProbePDGCode();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool ccprcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  bool ncprcok = proc_info.IsWeakNC() && (isP||isN) && (isnu||isnub);
  bool prcok   = ccprcok || ncprcok;
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
bool QELXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double Ethr = utils::kinematics::EnergyThreshold(interaction);
  if(E <= Ethr) {
     LOG("QELXSec", pINFO) << "Ev = " << E << " <= Ethreshold = "<< Ethr;
     return false;
  }
  return true;
}
//____________________________________________________________________________
void QELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::LoadConfig(void)
{
  //----- Get an algorithm to calculate differential cross sections dxsec/dQ2
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //----- Get an integrator
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *>(algf->GetAlgorithm(intgr));

  assert(fIntegrator);

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified
  //      (must be odd number for the Simpson rule)
  fNBins = fConfig->GetIntDef("N-logQ2-bins", 101);
  LOG("QELXSec", pDEBUG) << "Number of integration (logQ2) bins = " << fNBins;
  assert(fNBins>2);
}
//____________________________________________________________________________

