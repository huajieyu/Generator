#ifndef _COHRho_cohsigxy_H_
#define _COHRho_cohsigxy_H_

namespace genie  {
namespace utils {
namespace cohrhosigxy{
  //-- A realistic Breit-Wigner distribution with L-dependent width.
  double cohrhosigxy_cal(double ENU, double UU, double ELEP, double QQ2, double TMIN, double TMAX); 
} // cohrhosigxy namespace
} //utils namespace
} // genie  namespace
  
#endif   // _GENFUN_UTILS_H_
  
