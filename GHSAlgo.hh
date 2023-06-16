/// Implementation of the flavour dressing algorithm
/// by Rhorry Gauld, Alexander Huss and Giovanni Stagnitto,
/// as described in https://arxiv.org/abs/2208.11138 (v1)
///
/// Authors: Fabrizio Caola, Radoslaw Grabarczyk, Maxwell Hutt,
///          Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler

#ifndef __GHSALGO_HH__
#define __GHSALGO_HH__

#include "FlavInfo.hh"
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/NNH.hh"
#include "fastjet/PseudoJet.hh"
#endif

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

/// given a list of fiducial jets in the event (jets_in), return the
/// jets with "dressed" flavour information
///
/// @param jets_in: the input fiducial jets; the full list of event
/// particles and associated flavours will be deduced from
/// ClusterSequence associated with these jets. Note that a hardness cut
/// must have been applied to the set of input jets, otherwise the
/// algorithm is trivially IR unsafe.
///
/// @param beta: the value of beta used in the bottom-up soft-drop-like
/// phase for creating flavour clusters 
///
/// @param zcut: the value of zcut in the bottom-up soft-drop phase
///
/// @param Rcut: maximum 
///
/// @param alpha: power of (ktmax/ktmin) used in the flavour-kt distance
///
/// @returns the list of jets (in the same order as jets_in) dressed
/// with their flavour
std::vector<PseudoJet> run_GHS(const std::vector<PseudoJet> & jets_in, double ptcut = 20.0,
                            double beta=1, double zcut=0.1, double Rcut=0.1, double alpha=2,
                            double omega=0.0, bool use_fix_as2 = false);

std::vector<PseudoJet> create_flavoured_clusters_GHS(const std::vector<PseudoJet> & jets_in,
                                                     double beta=2, double zcut=0.1, double Rcut=0.1);

std::vector<PseudoJet> dress_GHS(const std::vector<PseudoJet> & jets_in, 
                                const std::vector<PseudoJet> & flavoured_clusters,
                                 double alpha, double omega, bool use_fix_as2 = false);


} // namespace contrib
FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif  // __GHSALGO_HH__
