#include "GHSAlgo.hh"
#include <limits>
// to facilitate use with fjcore
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "PseudoJetIO.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/NNH.hh"
#endif
#include "FlavInfo.hh"
#include <iostream>    
#include <iomanip>

//#define VERBOSE
//#define VVERBOSE

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

using namespace std;

const double _deltaR2_handover =
    pow(std::numeric_limits<double>::epsilon(), 0.5);

#ifdef VERBOSE
void print_PJ(ostream *ostr, const PseudoJet &p, unsigned precision,
              bool short_version, bool final_flav) {
  (*ostr).precision(precision);
  (*ostr) << setw(precision + 8) << p.px() << setw(precision + 8) << p.py()
          << setw(precision + 8) << p.pz() << setw(precision + 8) << p.E()
          << setw(precision + 8) << p.pt() << setw(precision + 8) << p.rap()
          << setw(precision + 8) << p.phi() << setw(precision + 8) << p.m()
          << setw(6) << p.user_index();
#ifdef PRINTMASSFLAVINFO          
  (*ostr) << setw(6) << p.has_user_info<MassFlavHistory>();
#else 
  (*ostr) << setw(6) << "";
#endif
  if (p.has_user_info<FlavHistory>()) {
    if (final_flav) (*ostr) << setw(14) << FlavHistory::current_flavour_of(p).description();
    else            (*ostr) << setw(14) << FlavHistory::initial_flavour_of(p).description();
    //cout << json(FlavHistory::current_flavour_of(p)) << endl;
  }
  if (!short_version) {
    (*ostr) << setw(6) << p.cluster_hist_index();
  }
}
#endif

// to use nnh when accumulating radiation to form flavoured clusters, with precise angles
class CamPrecBriefJet {
  public:
    void init(const PseudoJet & jet) {

      _pt2 = jet.pt2();
      double pt = sqrt(_pt2);
      _nx = jet.px() / pt;
      _ny = jet.py() / pt;
      _rap = jet.rap();
      _phi = jet.phi();
      constexpr double rap_transition = 0.1;
      if (fabs(_rap) < rap_transition) {
      // for moderately small rapidities switch to special rapidity formula
      //
      // rap = 1/2 log[(E+pz)/(E-pz)]
      //     = 1/2 log[1 + 2pz/(E-pz)]
      // 
      // and use log1p for the latter
        _rap = 0.5 * log1p(2*jet.pz()/(jet.E() - jet.pz()));
  
      }
   
    }
 
   double distance(const CamPrecBriefJet * other) const {
    double dphi = fabs(_phi - other->_phi);
    if (dphi > pi) {dphi = twopi - dphi;}
    double drap = _rap - other->_rap;

    // add special handling of dphi at small angles
    // (limits rounding errors when handling particles
    // close to the x or y axis)
    constexpr double phi_transition = 0.1;
    if (dphi < phi_transition) {
      // take a cross product of the n's (normalised), which 
      // is simply equal to sin(delta_phi)
      double cross = _nx * other->_ny - other->_nx * _ny;
      assert(cross >= -1.0 && cross <= 1.0);    
      // the sign can come out negative, but this isn't an issue
      // because we will use it in a context where the sign 
      // disappears
      dphi = asin(cross);
    }

    return sqrt(dphi*dphi + drap*drap);

   }
 
   double beam_distance() const {
     return numeric_limits<double>::max();
   }
  private:
    double _rap, _phi, _pt2, _nx, _ny;


};

struct GHSInfo {
  double alpha = 2.0;
  double omega = 0.0;
  unsigned njets;
  vector<PseudoJet> jets;
};

class GHSBriefJet {

public:
  void init(const PseudoJet & jet, GHSInfo * info) {
    _jet = jet;
    _info = info;
    _pt2 = jet.pt2();
    double pt = sqrt(_pt2);
    _nx = jet.px() / pt;
    _ny = jet.py() / pt;
    _rap = jet.rap();
    _phi = jet.phi();
    constexpr double rap_transition = 0.1;
    if (fabs(_rap) < rap_transition) {
      // for moderately small rapidities switch to special rapidity formula
      //
      // rap = 1/2 log[(E+pz)/(E-pz)]
      //     = 1/2 log[1 + 2pz/(E-pz)]
      // 
      // and use log1p for the latter
        _rap = 0.5 * log1p(2*jet.pz()/(jet.E() - jet.pz()));
  
    }

    // if it's a particle, see if it's associated with any of the jets
    // and record the information
    if (is_particle()) {
      for (unsigned i = 0; i < _info->jets.size(); i++) {
        //cout << jet << " associated with ?" << _info->jets[i] << ": ";
        if (jet.is_inside(_info->jets[i])) {
          _associated_jet = i;
          //cout << "yes" << endl;
          break;
        }
        //cout << "no" << endl;
      }
    }

    // now, evalute the beam distance for the PseudoJet
   if (is_jet()) {
     // no beam distance
     _diB = numeric_limits<double>::max();
   } else if (is_particle()) {

     // beam distance is defined only if the particle has no assigned jet
     if (_associated_jet == -1) {
       double ptBp = 0, ptBm = 0;
       for (const auto & j: _info->jets) {
         // _info->jets are PJs not BriefJets (as they must be to have an associated cs)
         // so we must recalculate the precise rapidity here.
         constexpr double rap_transition = 0.1;
         double jrap = j.rap(); if (fabs(jrap) < rap_transition) {jrap = 0.5 * log1p(2*j.pz()/(j.E() - j.pz()));}
         double jetrap = jet.rap(); if (fabs(jetrap) < rap_transition) {jetrap = 0.5 * log1p(2*jet.pz()/(jet.E() - jet.pz()));}
             double dyj = jrap - jetrap;
         ptBp += j.pt() * exp(min(0.0, -dyj));
         ptBm += j.pt() * exp(min(0.0, dyj));
       }
       double alpha = _info->alpha;
       double diBp =  max(pow(jet.pt(),   alpha), pow(ptBp,   alpha)) 
                     *min(pow(jet.pt(), 2-alpha), pow(ptBp, 2-alpha));
       double diBm =  max(pow(jet.pt(),   alpha), pow(ptBm,   alpha)) 
                     *min(pow(jet.pt(), 2-alpha), pow(ptBm, 2-alpha));
       _diB = min(diBp, diBm);
     } else {
       _diB = numeric_limits<double>::max();
     }

   } else {
     assert(false && "the PseudoJet should either be a particle or a jet, but was neither");

   }
  }
  bool   is_jet()        const {return _jet.user_index() == 1;}
  bool   is_particle()   const {return _jet.user_index() == 0;}

  static double geometrical_distance_squared(const GHSBriefJet *first, const GHSBriefJet *other){

    // do straight rapidity difference, because we already took
    // care of making the rapidity as accurate as possible
    double delta_y = first->_rap - other->_rap;
    double delta_phi = std::fabs(first->_phi - other->_phi);
    if (delta_phi > pi) delta_phi = twopi - delta_phi;

    // transition is somewhat arbitrary, but should be such that
    // we are in a region where arcsin() is unambiguous; can be
    // O(1), but must be < pi/2
    constexpr double phi_transition = 0.1;

    if (delta_phi < phi_transition) {
      // take a cross product of the n's (normalised), which 
      // is simply equal to sin(delta_phi)
      double cross = first->_nx * other->_ny - other->_nx * first->_ny;      
      // the sign can come out negative, but this isn't an issue
      // because we will use it in a context where the sign 
      // disappears
      delta_phi = asin(cross);
    }

    // omega = 0: we use deltaR^2 explicitly
    double omega = first->_info->omega;
    if (omega == 0.0) {
      return delta_y*delta_y + delta_phi*delta_phi;
    } else {
      double deltaR2 = delta_y*delta_y + delta_phi*delta_phi;
      if (deltaR2 > _deltaR2_handover)
        return 2*((cosh(omega*delta_y) - 1)/(omega*omega) - (cos(delta_phi) - 1));
      else
        return deltaR2;
    }
  }

  double distance(const GHSBriefJet * other_in) {
    // make sure that if at least one of them is a particle, then
    // the first is always a particle
    const GHSBriefJet * first = this;
    const GHSBriefJet * other = other_in;
    if (other->is_particle()) std::swap(first, other);

    // this can only happen if both are jets
    if (first->is_jet()) return numeric_limits<double>::max();

    // return the particle-particle distance (REVISIT FLAVOUR SUM?)
    if (other->is_particle()) {
      // give a flav dependence (only opposite flavours can annihilate)
      FlavInfo flavA = this->_jet.user_info<FlavHistory>().current_flavour();
      FlavInfo flavB = other_in->_jet.user_info<FlavHistory>().current_flavour();
      // if(!(flavA + flavB).is_flavourless()) {
      //   return numeric_limits<double>::max();
      // } else {
        return dij(first,other);
      // }
    } else {
      // ---- other must be a jet -----
      // first check to see if this particle is associated with any jet at all
      // (if not there is no distance)
      if (first->_associated_jet < 0) return numeric_limits<double>::max();
      // then see if it's associated with specifically with other 
      if (_info->jets[first->_associated_jet].cluster_hist_index() 
                                == other->_jet.cluster_hist_index()) {
        return dij(first,other);
      } else {
        return numeric_limits<double>::max();
      }
    }
  }
  double beam_distance() const {return _diB;}

  static double dij(const GHSBriefJet * first, const GHSBriefJet * other) {

    double alpha = first->_info->alpha;
    double ptf = sqrt(first->_pt2);
    double pto = sqrt(other->_pt2);
#ifdef VVERBOSE
        cout << " i (pt rap phi ui) = " << sqrt(first->_pt2) << " " << first->_rap << " " << first->_phi << " " << first->is_jet() << endl
             << " j (pt rap phi ui) = " << sqrt(other->_pt2) << " " << other->_rap << " " << other->_phi << " " << other->is_jet() << endl
             << " with distance = " << geometrical_distance_squared(first,other) 
            * max(pow(ptf,  alpha), pow(pto,  alpha))
            * min(pow(ptf,2-alpha), pow(pto,2-alpha))
             << endl;
#endif
    return geometrical_distance_squared(first,other) 
            * max(pow(ptf,  alpha), pow(pto,  alpha))
            * min(pow(ptf,2-alpha), pow(pto,2-alpha));

  }
  int associated_jet(){
    return _associated_jet;
  }
  void associate_with(int index_in_jets_array){
    _associated_jet = index_in_jets_array;
  }

private:
  int       _associated_jet = -1;
  PseudoJet _jet;
  GHSInfo * _info;
  double  _diB, _pt2, _rap, _phi, _nx, _ny;
};

// given initial state particles, it outputs the f^hat "aggregated flavour" jets as defined in the paper
// using the softdrop-like step, together with the information needed to associate them with a given jet.
std::vector<PseudoJet> create_flavoured_clusters_GHS(const std::vector<PseudoJet> &jets_in, double beta,
                                                   double zcut, double Rcut) {

  assert(jets_in.size() != 0);
  const ClusterSequence & cs = *(jets_in[0].associated_cs());
  vector<PseudoJet> initial_particles(cs.jets().begin(), cs.jets().begin() + cs.n_particles()); 
  int n_particles = initial_particles.size();
  NNH<CamPrecBriefJet> nnh(initial_particles);

  // if the user index has this value, then this particle has expired
  constexpr int expired_user_index = -100;
  // if it has this value, the particle is still relevant
  constexpr int current_user_index =  0;
  // make sure all initial particles are considered "current"
  for (auto & p: initial_particles) p.set_user_index(current_user_index);

  int iA, iB;
  while(n_particles > 1){

    double dij = nnh.dij_min(iA,iB);

    if (dij > Rcut) {
#ifdef VERBOSE
      cout << "dij > Rcut: break." << endl;
#endif
      break; // step 2 in the algorithm requires dij < Rcut.
    }

    FlavInfo flavA = initial_particles[iA].user_info<FlavHistory>().current_flavour();
    FlavInfo flavB = initial_particles[iB].user_info<FlavHistory>().current_flavour();
    // why do we need to call this?
    flavA.update_flavourless_attribute();
    flavB.update_flavourless_attribute();

#ifdef VERBOSE
    cout << "dij = " << dij << ", iA = " << iA << ", iB = " << iB << endl;
    print_PJ(&cout, initial_particles[iA], 5, true, true); cout << endl;
    print_PJ(&cout, initial_particles[iB], 5, true, true); cout << endl;
#endif

    // attributes that we will use below
    // NB: default for i_flavoured and i_flavourless
    //     are used even when the particles are not flavoured
    int i_flavoured   = iA; 
    int i_flavourless = iB; 
    bool both_flavoured = false;
    bool both_flavourless = false;

    if (!flavA.is_flavourless() && !flavB.is_flavourless()) {
      both_flavoured = true;
    } else if (!flavB.is_flavourless()) {
      //do_merge(iA,iB);
      // together with first failure of first if condition, this implies
      // that flavA is flavourless
      i_flavoured   = iB;
      i_flavourless = iA;
    } else if (flavA.is_flavourless()) { 
      //flavB is flavourless, so if flavA is flavourless, both are flavourless
      both_flavourless = true;
    } else {
      // if none of these hold, then flavA is flavoured, and flavB is
      // flavourless; both_flavoured and both_flavourless are false.
      // Our default iA and iB settings above are good
    }

    // depending on the flavours, we perform, or do not perform softdrop:
    if(both_flavoured) {
#ifdef VERBOSE
      cout << "both flavoured. accumulation for these is complete." << endl;
#endif
      // the fix that was communicated in les Houches is that
      // in case both pseudojets are flavoured, we add them as well (with net
      // flavour summation)
      // NB: not sure this works!
      // if (use_fix_as2) {
      //   // check SD condition
      //   double z = min(initial_particles[iA].pt(), initial_particles[iB].pt())
      //               /(initial_particles[iA].pt() + initial_particles[iB].pt());

      //   // there is an issue with pow() in the dddreal/qdreal numerical types
      //   // The types use pow(a,b) = e^(b*log(a)).
      //   // Instead we differentiate cases as follows:
      //   //double cut = zcut*pow(dij/Rcut, beta);
      //   assert(beta > 0);
      //   double cut = zcut*pow(dij/Rcut, beta);

      //   // if SD condition passes, we merge the two flavoured particles
      //   if (z > cut) {
      //     // merge
      //     PseudoJet new_pseudojet = initial_particles[iA];
      //     new_pseudojet.reset_momentum(initial_particles[iA] + initial_particles[iB]); //<- resetting only the momentum keeps the
      //     new_pseudojet.set_user_index(current_user_index);
      //     FlavInfo flav = FlavHistory::current_flavour_of(initial_particles[iA]) + FlavHistory::current_flavour_of(initial_particles[iB]);
      //     /// set FlavInfo attribute
      //     new_pseudojet.set_user_info(new FlavHistory(flav));

      //     initial_particles.push_back(new_pseudojet);
      //     nnh.merge_jets(iA, iB, new_pseudojet, initial_particles.size()-1);
      //     initial_particles[iA].set_user_index(expired_user_index);
      //     initial_particles[iB].set_user_index(expired_user_index);
      //   } else {
      //     // Otherwise, remove them like in the original algorithm
      //     nnh.remove_jet(iA);
      //     nnh.remove_jet(iB);
      //     n_particles -= 2;
      //     continue;
      //   }
      // } else { // < original algorithm as in 2208.11138v1
        nnh.remove_jet(iA);
        nnh.remove_jet(iB);
        n_particles -= 2;
        continue;
      // }

    } else {

      double z = min(initial_particles[iA].pt(), initial_particles[iB].pt())
                   /(initial_particles[iA].pt() + initial_particles[iB].pt());

      // there is an issue with pow() in the dddreal/qdreal numerical types
      // The types use pow(a,b) = e^(b*log(a)).
      // Instead we differentiate cases as follows:
      //double cut = zcut*pow(dij/Rcut, beta);
      assert(beta > 0);
      double cut = zcut*pow(dij/Rcut, beta);

      if ( (z > cut) || both_flavourless ){
#ifdef VERBOSE
      cout << "both flavourless, or SD cut. merge iA, iB" << endl;
#endif
        // if softdrop condition is passed, or both particles are flavourless, merge
        PseudoJet new_pseudojet = initial_particles[i_flavoured];
        new_pseudojet.reset_momentum(initial_particles[iA] + initial_particles[iB]); //<- resetting only the momentum keeps the
        new_pseudojet.set_user_index(current_user_index);
        // information about the jet to which the flavour is associated.
        //initial_particles.erase(initial_particles.begin() + i_flavourless);
        //initial_particles.erase(initial_particles.begin() + i_flavoured);
        initial_particles.push_back(new_pseudojet);
        nnh.merge_jets(iA, iB, new_pseudojet, initial_particles.size()-1);
        initial_particles[i_flavoured  ].set_user_index(expired_user_index);
        initial_particles[i_flavourless].set_user_index(expired_user_index);

      } else {
#ifdef VERBOSE
      cout << "SD criterion doesn't pass. Removing that flavourless object" << endl;
#endif
        // when one is flavoured, and the softdrop condition is not
        // passed, we remove the flavourless particle
        nnh.remove_jet(i_flavourless);
        //initial_particles.erase(initial_particles.begin() + i_flavourless);
        initial_particles[i_flavourless].set_user_index(expired_user_index);

      }
    }
    n_particles--;

  } // end of while loop
  vector<PseudoJet> flavoured_clusters;
  for (auto & a : initial_particles) {
    const FlavInfo & flav = a.user_info<FlavHistory>().current_flavour();
    if (a.user_index() != expired_user_index && (!(flav.is_flavourless())) ){
      flavoured_clusters.push_back(a);
#ifdef VERBOSE
    cout << "accumulated flavour jet : "; print_PJ(&cout, a, 8, true, true); cout << " " << a.has_valid_cluster_sequence() << endl;
#endif
    }
  }

  return flavoured_clusters;

}

std::vector<PseudoJet> dress_GHS(const std::vector<PseudoJet> & jets_in, 
                                const std::vector<PseudoJet> & flavoured_clusters,
                                 double alpha, double omega, bool use_fix_as2){
  vector<PseudoJet> clusters = flavoured_clusters;
  vector<PseudoJet> all;
  vector<PseudoJet> jets = jets_in;
  vector<FlavInfo> jet_flavs(jets.size(), 0);
  int njets = jets.size();

  GHSInfo ghs_info;
  ghs_info.jets  = jets;
  ghs_info.njets = jets.size();
  ghs_info.alpha = alpha;
  ghs_info.omega = omega;
  
  // make sure that jets have no flavour, and give them an index meaning that they are a jet
  for (auto & j : jets) {
    j.set_user_info(new FlavHistory(FlavInfo( 0)));
    //FlavInfo flavj = j.user_info<FlavHistory>().current_flavour();
    //flavj.update_flavourless_attribute();
    j.set_user_index(1); //< the signal that it is a jet
    all.push_back(j);
  }
  for (auto & c : clusters){
    c.set_user_index(0); //< the signal that it is a flavour cluster
    all.push_back(c);
  }
  

  if (all.size() == 0) {return all;}

#ifdef VERBOSE
  cout << "initial jets + clusters: (0 = cluster, 1 = jet)" << endl;
  for (auto & a : all){
    print_PJ(&cout, a, 5, true, true); cout << endl;
  }
#endif
  //all.insert(all.end(), clusters.begin(), clusters.end());
  // set up nnh
  NNH<GHSBriefJet,GHSInfo> nnh(all, &ghs_info);
  int iA, iB;
  while (njets > 0) { //the loop does not change njets, but njets > 0 is necessary given that the selector could cut off all jets

    double dij = nnh.dij_min(iA,iB);
#ifdef VERBOSE
    cout << iA << " " << iB << endl;
    cout << "dij = "<< dij << endl;
#endif

    // LS-2023-02-10: not sure this is very safe...
    if (dij == numeric_limits<double>::max()) {
      break;
    }
    if (iB >= 0) {
      if (iA > iB) std::swap(iA,iB);
      // we must never have two jets
      assert (iB >= njets && "second entry must be a particle");
      if (iA < njets) {
        // if the first is a jet, assign B's flavour to A and then remove B
        // (note that through the shared pointer, this also affects the
        // flavour of the objects in the NNH object -- which is dangerous -- one
        // should really remove the jet and add it back in)
        //FlavInfo * flav_info = dynamic_cast<FlavInfo*>(jets[iA].user_info_shared_ptr().get());
        //*flav_info = *flav_info + all[iB].user_info<FlavInfo>();
        //FlavInfo flavA = jets[iA].user_info<FlavHistory>().current_flavour();
        FlavInfo flavB = all[iB].user_info<FlavHistory>().current_flavour();
        jet_flavs[iA] = jet_flavs[iA] + flavB;
        //jets[iA].set_user_info(new FlavHistory(flavA + flavB));
        nnh.remove_jet(iB);
#ifdef VERBOSE
        cout << "adding flavour of " << iB << " to " << iA << ", removing " << iB << endl;
#endif
      } else {
        if (use_fix_as2) {
          // merge
          PseudoJet new_pseudojet = all[iA];
          new_pseudojet.reset_momentum(all[iA] + all[iB]); //<- resetting only the momentum keeps the
          new_pseudojet.set_user_index(0);
          FlavInfo flav = FlavHistory::current_flavour_of(all[iA]) + FlavHistory::current_flavour_of(all[iB]);
          /// set FlavInfo attribute
          new_pseudojet.set_user_info(new FlavHistory(flav));
#ifdef VERBOSE
        cout << "two flavoured clusters combine." << endl;
        cout << "iA = "  << iA; print_PJ(&cout, all[iA], 12, true, true); cout << endl;
        cout << "iB = "  << iB; print_PJ(&cout, all[iB], 12, true, true); cout << endl;
        cout << "into "; print_PJ(&cout, new_pseudojet, 12, true, true); cout << endl;
#endif          
          all.push_back(new_pseudojet);
          nnh.merge_jets(iA, iB, new_pseudojet, all.size()-1);
        } else {
        // both are particles, so annihilate their flavour, i.e. remove them from the list
#ifdef VERBOSE
        cout << "two flavoured clusters annihilate." << endl;
#endif
        nnh.remove_jet(iA);
        nnh.remove_jet(iB);
        }
      }
    } else {
   //   assert(iA >= njets && "for beam clustering, iA must be a particle");
#ifdef VERBOSE
      cout << "removing " << iA << endl;
#endif
      nnh.remove_jet(iA);
    }

  }
  for (unsigned i = 0; i < jets.size(); i++){
    jet_flavs[i].update_flavourless_attribute();
    jets[i].set_user_info(new FlavHistory(FlavInfo(jet_flavs[i])));
    // restore user index to what it was 
    jets[i].set_user_index(jets_in[i].user_index());
  }
  return jets;
}


LimitedWarning ghs_fiducial_warning(10);

//---------------------------------------------------------------------------
// this contains the flavour clusters + the dressing. in a separate class to have space for other association criterions (to add)
std::vector<PseudoJet> run_GHS(const std::vector<PseudoJet> & jets_in, double ptcut,
                              double beta, double zcut, double Rcut, double alpha, double omega,
                              bool use_fix_as2) {

  if (jets_in.size() == 0) return jets_in;

#ifdef VERBOSE
  cout << " -- generate flavour clusters -- " << endl;
#endif
  /// these are the "flavoured particles and clusters" from 2208.11138v1, p.2
  vector<PseudoJet> flavoured_clusters = create_flavoured_clusters_GHS(jets_in, beta, zcut, Rcut);

  /// different (non-default) association criteria would go here. For instance (in the future): cluster flavoured_clusters as ghosts
  /// with jets_in; or test whether a cluster is less than R_tag away from a jet. By default, flavoured_clusters are associated
  /// as said in the paper, by the association of the original flavoured particle to an anti-kt jet during the earlier flavour-agnostic
  /// clustering.

  ghs_fiducial_warning.warn("Watch out: in run_GHS a fiducial selector is being applied to jets_in (this should be removed; but then run_GHS needs to get a list of all particles)");
  Selector select_pt = SelectorPtMin(ptcut); 
  vector<PseudoJet> selected_jets = select_pt(jets_in);

  if (selected_jets.size() == 0) {
    return selected_jets;
  }

#ifdef VERBOSE
  cout << " -- dressing jets -- " << endl;
#endif
  vector<PseudoJet> final_jets = dress_GHS(selected_jets, flavoured_clusters, alpha, omega, use_fix_as2);
#ifdef VERBOSE
  cout << "final dressed jets = " << endl;
  for (auto & a : final_jets){
    print_PJ(&cout, a, 5, true, true); cout << endl;
  }
#endif
  
  return final_jets;
}

} // namespace contrib
FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
