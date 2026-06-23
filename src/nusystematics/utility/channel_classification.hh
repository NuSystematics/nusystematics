#pragma once

// Per-event channel classification used for dial-applicability filtering.
//
// A channel key is a string of the form
//
//     <CC|NC>_<nu-species>_<target-name>_<topology>
//
//     e.g.  "CC_numu_Ar40_QE"
//           "NC_nuebar_C12_RES"
//           "CC_numu_Ar40_Other"     (unknown topology bucket)
//
// `applies_to_channels` is a vector<string> of fnmatch-style glob patterns
// matched against the channel key. Examples:
//
//   ["CC_*_*_QE"]       -> all CC quasi-elastic events
//   ["CC_numu_*"]       -> all numu CC events
//   ["*_QE"]            -> QE on any current / nu / target
//   []                  -> applies to everything (no filtering)
//
// The pattern matcher is POSIX fnmatch with default flags, so `*` matches any
// sequence (including underscores) and `?` matches a single character.

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

#include <cstdio>
#include <cstdlib>
#include <fnmatch.h>
#include <string>
#include <vector>

namespace nusyst::channel {

inline std::string NuName(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "numu"  : "numubar";
    case 12: return pdg > 0 ? "nue"   : "nuebar";
    case 16: return pdg > 0 ? "nutau" : "nutaubar";
    default: {
      char buf[16];
      std::snprintf(buf, sizeof(buf), "nu%d", pdg);
      return buf;
    }
  }
}

inline std::string TargetName(int A) {
  switch (A) {
    case 40: return "Ar40";
    case 12: return "C12";
    case 16: return "O16";
    case 56: return "Fe56";
    case 1:  return "H";
    default: {
      char buf[16];
      std::snprintf(buf, sizeof(buf), "A%d", A);
      return buf;
    }
  }
}

inline std::string Topology(const genie::EventRecord &ev) {
  auto const &proc = ev.Summary()->ProcInfo();
  if (proc.IsQuasiElastic())       return "QE";
  if (proc.IsMEC())                return "MEC";
  if (proc.IsResonant())           return "RES";
  if (proc.IsDeepInelastic())      return "DIS";
  if (proc.IsCoherentProduction()) return "COH";
  if (proc.IsCoherentElastic())    return "COH";
  return "Other";
}

inline std::string Current(const genie::EventRecord &ev) {
  return ev.Summary()->ProcInfo().IsWeakCC() ? "CC" : "NC";
}

inline std::string MakeChannelKey(const genie::EventRecord &ev) {
  int nu_pdg = ev.Probe() ? ev.Probe()->Pdg() : 0;
  int tgt_A  = ev.TargetNucleus() ? ev.TargetNucleus()->A() : 1;
  return Current(ev) + "_" + NuName(nu_pdg) + "_" + TargetName(tgt_A) + "_" + Topology(ev);
}

// Returns true if `channel_key` matches at least one fnmatch glob in
// `patterns`. An empty pattern list is treated as "match everything"
// (un-filtered provider).
inline bool MatchesAny(const std::string &channel_key,
                       const std::vector<std::string> &patterns) {
  if (patterns.empty()) return true;
  for (auto const &p : patterns) {
    if (::fnmatch(p.c_str(), channel_key.c_str(), 0) == 0) return true;
  }
  return false;
}

} // namespace nusyst::channel
