#ifndef __PSEUDOJETIO_HH__
#define __PSEUDOJETIO_HH__
#include <iostream>
#include "fastjet/PseudoJet.hh"
#include "FlavInfo.hh"

inline std::ostream & operator<<(std::ostream & ostr, const fastjet::PseudoJet & p) {
  ostr <<   "pt=" << p.pt();
  ostr << ", y=" << p.rap();
  ostr << ", phi=" << p.phi();
  ostr << ", m=" << p.m();
  if (p.has_user_info<fastjet::contrib::FlavInfo>()) {
      ostr << ", " << p.user_info<fastjet::contrib::FlavInfo>().description();
  }
  return ostr;
}


#endif // __PSEUDOJETIO_HH__
