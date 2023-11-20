/// Implementation of the flavour dressing algorithm
/// by Rhorry Gauld, Alexander Huss, and Giovanni Stagnitto,
/// as described in https://arxiv.org/abs/2208.11138 (v2)
///
/// initial implementation by:
///   Fabrizio Caola, Radoslaw Grabarczyk, Maxwell Hutt,
///   Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
/// adapted to v2 by:
///   Rhorry Gauld, Alexander Huss, and Giovanni Stagnitto


#ifndef __GHSALGO_HH__
#define __GHSALGO_HH__

#include "FlavInfo.hh"
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/NNH.hh"
#include "fastjet/PseudoJet.hh"
#endif

FASTJET_BEGIN_NAMESPACE  // defined in fastjet/internal/base.hh
namespace contrib {

  /// given a list of base-jets (before applying a hardness cut) in
  /// the event (jets_base), return the jets with "dressed" flavour information
  ///
  /// @param jets_base: the input base-jets; the full list of event
  /// particles and associated flavours will be deduced from
  /// ClusterSequence associated with these jets.
  ///
  /// @param ptcut: this parameter applies a hardness cut on the base-jets
  /// and should be chosen to match the fiducial jet definition
  ///
  /// @param alpha: power of (ktmax/ktmin) used in the flavour-kt distance
  ///
  /// @param omega: relative weighting of rapidity separation
  ///
  /// @returns the list of hard jets dressed with their flavour
  std::vector<PseudoJet> run_GHS(
      const std::vector<PseudoJet> &jets_base, double ptcut, double alpha = 1.,
      double omega = 2.,
      const FlavRecombiner &flav_recombiner = FlavRecombiner());

}  // namespace contrib
FASTJET_END_NAMESPACE  // defined in fastjet/internal/base.hh

#endif  // __GHSALGO_HH__
