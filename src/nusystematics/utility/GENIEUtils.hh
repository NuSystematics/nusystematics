#pragma once

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/utility/simbUtility.hh"

#include "Framework/EventGen/EventRecord.h"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"

#include "Framework/Interaction/SppChannel.h"
#include "Framework/Interaction/ProcessInfo.h"

#include <sstream>

namespace nusyst {
/// Gets the GENIE SPP channel enum for a supplied GHepEvent
///
/// N.B. going via NEUT mode as genie::SppChannel::FromInteraction doesn't work
/// for RES events.
inline genie::SppChannel_t SPPChannelFromGHep(genie::EventRecord const &ev) {
  int NEUTCh = genie::utils::ghep::NeutReactionCode(&ev);

  switch (NEUTCh) {
  case 11: { // kCCSPP_PPip
    return genie::kSpp_vp_cc_10100;
  }
  case 12: { // kCCSPP_PPi0
    return genie::kSpp_vn_cc_10010;
  }
  case 13: { // kCCSPP_NPip
    return genie::kSpp_vn_cc_01100;
  }
  case -11: {
    return genie::kSpp_vbn_cc_01001;
  }
  case -12: {
    return genie::kSpp_vbp_cc_01010;
  }
  case -13: {
    return genie::kSpp_vbp_cc_10001;
  }
  default: {
    return genie::kSppNull;
  }
  }
}

enum class QELikeTarget_t { kNN = 0, knp, kQE, kInvalidTopology };
NEW_SYSTTOOLS_EXCEPT(indeterminable_QELikeTarget);


// TH: Adapted IsPrimary function from NUISANCE
inline bool IsPrimary(genie::EventRecord const &ev, genie::GHepParticle *p){
  
  if (p->Status() == genie::kIStInitialState || p->Status() == genie::kIStNucleonTarget){
        return true;
      }

      // Reject intermediate states
      if (p->Status() == genie::kIStDISPreFragmHadronicState || 
        p->Status() == genie::kIStPreDecayResonantState || 
        p->Status() == genie::kIStInitialState || 
        p->Status() == genie::kIStDecayedState || 
        p->Status() == genie::kIStUndefined){ 
          return false;
      }
      
      // Check if the mother is the neutrino or IS nucleon
      if (p->FirstMother() < 2){
        return true;
      }

      // Loop over particle's mothers, gmothers... and clean out intermediate particles
      genie::GHepParticle *mother = ev.Particle(p->FirstMother());
      while (mother->FirstMother() > 1) {
        // could be mother's status is actually a decayed state linked back to the vertex
        if (mother->Status() == genie::kIStDecayedState || // a decayed state
            mother->Status() == genie::kIStDISPreFragmHadronicState ||  // a DIS state before fragementation
            mother->Status() == genie::kIStPreDecayResonantState) { // a pre-decay resonant state
              mother = ev.Particle(mother->FirstMother());
            } else { // if not move out of the loop
              break;
            }
      }

      // Then check is mother is associated with primary
      int MotherID = mother->FirstMother();
      if (MotherID > 2){
        return false;
      }

      // Finally, this could mean particle is marked for transport 
      // Could also be interactions of a free proton
      if (p->Status() == genie::kIStHadronInTheNucleus ||  // Then require the particle to be paseed to FSI
      (p->Status() == genie::kIStStableFinalState && // Can also have interaction on free proton
       ev.Summary()->InitState().TgtPtr()->A() == 1 &&
       ev.Summary()->InitState().TgtPtr()->Z() == 1) ) {
          return true;
      } 
  
    return false;
  
}

inline QELikeTarget_t GetQELikeTarget(genie::EventRecord const &ev) {

  if (ev.Summary()->ProcInfo().IsQuasiElastic() &&
      !ev.Summary()->ExclTag().IsCharmEvent()) {
    return QELikeTarget_t::kQE;
  }
  if (ev.Summary()->ProcInfo().IsMEC()) {
    genie::Target const &tgt = ev.Summary()->InitState().Tgt();
    size_t nuc_pdg = tgt.HitNucPdg();

    switch (nuc_pdg - 2000000200) {
    case 0:
    case 2: {
      return QELikeTarget_t::kNN;
    }
    case 1: {
      return QELikeTarget_t::knp;
    }
    default: {
      throw indeterminable_QELikeTarget()
          << "[ERROR]: Failed to determine 2p2h topology from interaction: "
          << ev.Summary()->AsString()
          << ", expected the target nucleon PDG: " << nuc_pdg
          << " - 2000000200 = " << (nuc_pdg - 2000000200)
          << " to be 0, 1, or 2.";
    }
    }
  }

  return QELikeTarget_t::kInvalidTopology;
}

/// Encoded as IsNeutrino * 1 + IsCC * 10 + TargetIsProton * 100 + NPi * 1000 +
/// NPiplus * 10000 + NPiminus * 100000 + NPi0 * 1000000
///
///\n.b. is TargetIsProton == 2 then it can be either proton or neutron target
typedef size_t NRPiChan_t;

inline bool IsNeutrinoNRPiChan(NRPiChan_t ch) { return ch % 10; }
inline bool IsCCNRPiChan(NRPiChan_t ch) { return (ch / 10) % 10; }
inline bool IsProtonTargetNRPiChan(NRPiChan_t ch) {
  return ((ch / 100) % 10) & 2;
}
inline bool IsNeutronTargetNRPiChan(NRPiChan_t ch) {
  return ((ch / 100) % 10) & 1;
}
inline size_t GetNRPiChanNPi(NRPiChan_t ch) { return (ch / 1000) % 10; }
inline size_t GetNRPiChanNPip(NRPiChan_t ch) { return (ch / 10000) % 10; }
inline size_t GetNRPiChanNPim(NRPiChan_t ch) { return (ch / 100000) % 10; }
inline size_t GetNRPiChanNPi0(NRPiChan_t ch) { return (ch / 1000000) % 10; }

inline NRPiChan_t BuildNRPiChannel(bool IsNeutrino, bool IsCC,
                                   int TargetNucleon, size_t NPi,
                                   size_t NPiplus = 0, size_t NPiminus = 0,
                                   size_t NPi0 = 0) {
  return IsNeutrino * 1 + IsCC * 10 + TargetNucleon * 100 + NPi * 1000 +
         NPiplus * 10000 + NPiminus * 100000 + NPi0 * 1000000;
}

inline NRPiChan_t GetNRPiChannel(genie::EventRecord const &ev) {
  // This code in this method is adapted from the GENIE source code found in
  // GHep/GHepUtils.cxx This method therefore carries the GENIE licence as
  // copied below:
  //
  /// Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
  /// For the full text of the license visit http://copyright.genie-mc.org
  /// or see $GENIE/LICENSE
  //
  if (!ev.Summary()->ProcInfo().IsDeepInelastic()) {
    return 0;
  }

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  if (!tgt.HitNucIsSet()) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to get hit nucleon kinematics as it was not "
           "included in this GHep event. This is a fatal error.";
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }
  int NPi0 = 0, NPip = 0, NPim = 0;

  bool nuclear_target = tgt.IsNucleus();

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    genie::GHepStatus_t ghep_ist = (genie::GHepStatus_t)p->Status();
    int ghep_pdgc = p->Pdg();
    int ghep_fm = p->FirstMother();
    int ghep_fmpdgc = (ghep_fm == -1) ? 0 : ev.Particle(ghep_fm)->Pdg();

    // For nuclear targets use hadrons marked as 'hadron in the nucleus'
    // which are the ones passed in the intranuclear rescattering
    // For free nucleon targets use particles marked as 'final state'
    // but make an exception for decayed pi0's,eta's (count them and not their
    // daughters)

    bool decayed =
        (ghep_ist == genie::kIStDecayedState &&
         (ghep_pdgc == genie::kPdgPi0 || ghep_pdgc == genie::kPdgEta));
    bool parent_included =
        (ghep_fmpdgc == genie::kPdgPi0 || ghep_fmpdgc == genie::kPdgEta);

    bool count_it =
        (nuclear_target && ghep_ist == genie::kIStHadronInTheNucleus) ||
        (!nuclear_target && decayed) ||
        (!nuclear_target && ghep_ist == genie::kIStStableFinalState &&
         !parent_included);

    if (!count_it) {
      continue;
    }
    if (ghep_pdgc == genie::kPdgPiP) {
      NPip++;
    } else if (ghep_pdgc == genie::kPdgPiM) {
      NPim++;
    } else if (ghep_pdgc == genie::kPdgPi0) {
      NPi0++;
    }
  }

  return BuildNRPiChannel(ISLep->Pdg() > 0, ev.Summary()->ProcInfo().IsWeakCC(),
                          (tgt.HitNucPdg() == genie::kPdgProton) ? 2 : 1,
                          NPi0 + NPip + NPim, NPip, NPim, NPi0);
}

inline std::string GetNRPiChannelName(NRPiChan_t ch) {
  std::stringstream ss("");

  ss << "NR_" << (IsNeutrinoNRPiChan(ch) ? "nu" : "nubar") << "_"
     << (IsNeutronTargetNRPiChan(ch) ? "n" : "")
     << (IsProtonTargetNRPiChan(ch) ? "p" : "") << "_"
     << (IsCCNRPiChan(ch) ? "CC" : "NC") << "_" << GetNRPiChanNPi(ch) << "Pi";

  return ss.str();
}

/// Check if an event channel is equivalent to an overall channel, n and p
/// target events are equivalent to np channels and NPions > NMax pions are
/// equivalent to NPions == NMax pions
inline bool ChannelsAreEquivalent(NRPiChan_t ch, NRPiChan_t event_ch,
                                  size_t NMaxPions) {
  if (IsNeutrinoNRPiChan(ch) != IsNeutrinoNRPiChan(event_ch)) {
    return false;
  }
  if (IsCCNRPiChan(ch) != IsCCNRPiChan(event_ch)) {
    return false;
  }
  if (!((IsProtonTargetNRPiChan(event_ch) && IsProtonTargetNRPiChan(ch)) ||
        (IsNeutronTargetNRPiChan(event_ch) && IsNeutronTargetNRPiChan(ch)))) {
    return false;
  }
  size_t NPiChannel = (GetNRPiChanNPi(event_ch) > NMaxPions)
                          ? NMaxPions
                          : GetNRPiChanNPi(event_ch);
  if (NPiChannel != GetNRPiChanNPi(ch)) {
    return false;
  }
  return true;
}

inline double GetErecoil_MINERvA_LowRecoil(genie::EventRecord const &ev) {
  // Get total energy of hadronic system.
  double Erecoil = 0.0;

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    if (p->Status() != genie::kIStStableFinalState) {
      continue;
    }
    switch (p->Pdg()) {
    case 2212:
    case 211:
    case -211: {
      Erecoil += p->KinE();
      break;
    }
    case 111:
    case 11:
    case -11:
    case -22: {
      Erecoil += p->E();
      break;
    }
    default: {
    }
    }
  }

  // For nue CC scattering, we would have counted the E of the charged lepton,
  // subtract it off here
  if (ev.Summary()->ProcInfo().IsWeakCC() && (abs(ev.Probe()->Pdg()) == 12)) {
    Erecoil -= ev.FinalStatePrimaryLepton()->P4()->E();
  }

  return Erecoil;
}

inline simb_mode_copy GetSimbMode(genie::EventRecord const &ev) {

  simb_mode_copy mode = simb_mode_copy::kUnknownInteraction;

  if (ev.Summary()->ProcInfo().IsQuasiElastic()) {
    mode = simb_mode_copy::kQE;
  } else if (ev.Summary()->ProcInfo().IsMEC()) {
    mode = simb_mode_copy::kMEC;
  } else if (ev.Summary()->ProcInfo().IsResonant()) {
    mode = simb_mode_copy::kRes;
  } else if (ev.Summary()->ProcInfo().IsDeepInelastic()) {
    mode = simb_mode_copy::kDIS;
  } else if (ev.Summary()->ProcInfo().IsCoherentProduction()) {
    mode = simb_mode_copy::kCoh;
  } else if (ev.Summary()->ProcInfo().IsCoherentElastic()) {
    mode = simb_mode_copy::kCohElastic;
  } else if (ev.Summary()->ProcInfo().IsElectronScattering()) {
    mode = simb_mode_copy::kElectronScattering;
  } else if (ev.Summary()->ProcInfo().IsIMDAnnihilation()) {
    mode = simb_mode_copy::kIMDAnnihilation;
  } else if (ev.Summary()->ProcInfo().IsInverseBetaDecay()) {
    mode = simb_mode_copy::kInverseBetaDecay;
  } else if (ev.Summary()->ProcInfo().IsGlashowResonance()) {
    mode = simb_mode_copy::kGlashowResonance;
  } else if (ev.Summary()->ProcInfo().IsAMNuGamma()) {
    mode = simb_mode_copy::kAMNuGamma;
  } else if (ev.Summary()->ProcInfo().IsDiffractive()) {
    mode = simb_mode_copy::kDiffractive;
  } else if (ev.Summary()->ProcInfo().IsEM()) {
    mode = simb_mode_copy::kEM;
  } else if (ev.Summary()->ProcInfo().IsWeakMix()) {
    mode = simb_mode_copy::kWeakMix;
  }

  return mode;
}

inline std::string DumpGENIEEv(genie::EventRecord const &ev) {
  std::stringstream ss("");
  ss << ev.Summary()->AsString() << std::endl;

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();

  bool nuclear_target = tgt.IsNucleus();

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;
  size_t p_it = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    genie::GHepStatus_t ghep_ist = (genie::GHepStatus_t)p->Status();
    int ghep_pdgc = p->Pdg();
    int ghep_fm = p->FirstMother();
    int ghep_fmpdgc = (ghep_fm == -1) ? 0 : ev.Particle(ghep_fm)->Pdg();

    bool decayed =
        (ghep_ist == genie::kIStDecayedState &&
         (ghep_pdgc == genie::kPdgPi0 || ghep_pdgc == genie::kPdgEta));
    bool parent_included =
        (ghep_fmpdgc == genie::kPdgPi0 || ghep_fmpdgc == genie::kPdgEta);

    bool count_it =
        (nuclear_target && ghep_ist == genie::kIStHadronInTheNucleus) ||
        (!nuclear_target && decayed) ||
        (!nuclear_target && ghep_ist == genie::kIStStableFinalState &&
         !parent_included);

    ss << "Part: " << p_it++ << ", pdg = " << ghep_pdgc << ", decayed ? "
       << decayed << ", parent_included ? " << parent_included
       << ", counting ? " << count_it << std::endl;
  }
  ss << std::endl;
  return ss.str();
}

// Copy of https://github.com/NuSoftHEP/nugen/blob/6bcd82d9310bd0480df9a8ed03bfc8d8a2c80eff/nugen/EventGeneratorBase/GENIE/GENIE2ART.cxx#L98-L126
inline std::string ExpandEnvVar(const std::string& s)
{
  // utility function:
  // if input "s" starts w/ $, return corresponding env var value,
  // otherwise return as is

  if ( s.find('$') != 0 ) return s;

  // need to remove ${}'s
  std::string sEnvVar = s;
  char rmchars[] = "$(){} ";
  for (unsigned int i = 0; i < strlen(rmchars); ++i) {
    // remove moves matching characters in [first,last) to end and
    //   returns a past-the-end iterator for the new end of the range [funky!]
    // erase actually trims the string
    sEnvVar.erase( std::remove(sEnvVar.begin(), sEnvVar.end(),
                               rmchars[i]), sEnvVar.end() );
  }
  const char* charEnvValue = std::getenv(sEnvVar.c_str());
  if ( ! charEnvValue ) {
    // resolved into an empty string ... not what one would expect
      std::cout << " can't resolve " << s << " via getenv(\"" << sEnvVar << "\")" << std::endl;
    return s; // return original (though we won't get here due to throw)
  }
  return std::string(charEnvValue);

}

} // namespace nusyst
