#pragma once

// Named event variables calculated from a genie::EventRecord.
//
// Maps a variable name given in a YAML configuration (e.g. a BDT input feature
// or an event selection cut) to the calculation on the GENIE event record.
// To support a new variable, add one entry to BuildEventVariableRegistry().
//
// Units: GENIE native, i.e. energies and momenta in GeV.
// PhiAdler is in degrees, as in NUISANCE.
// A variable that is not defined for the event (e.g. no such particle) returns
// kUndefined.

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"
#include "nusystematics/utility/exceptions.hh"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Interaction/Interaction.h"

#include "systematicstools/utility/string_parsers.hh"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace nusyst {
namespace evvar {

NEW_SYSTTOOLS_EXCEPT(unknown_event_variable);
NEW_SYSTTOOLS_EXCEPT(invalid_event_cut);

using EventVariable = std::function<double(genie::EventRecord const &)>;

constexpr double kUndefined = -999.;

// Nucleon mass used for W_nuc_rest and x_Bj; same value as NUISANCE
constexpr double kNucleonMass = 0.93827203; // GeV

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

inline TLorentzVector ProbeP4(genie::EventRecord const &ev) {
  return *ev.Probe()->P4();
}

inline bool HasFSLepton(genie::EventRecord const &ev) {
  return ev.FinalStatePrimaryLepton() != nullptr;
}

inline TLorentzVector FSLeptonP4(genie::EventRecord const &ev) {
  return *ev.FinalStatePrimaryLepton()->P4();
}

// Four-momentum transfer, q = p_nu - p_lep
inline TLorentzVector MomentumTransfer(genie::EventRecord const &ev) {
  return ProbeP4(ev) - FSLeptonP4(ev);
}

// Hadrons from the primary vertex, before FSI
inline bool IsPreFSIHadron(genie::EventRecord const &ev,
                           genie::GHepParticle *p) {
  if (p->Status() == genie::kIStInitialState ||
      p->Status() == genie::kIStNucleonTarget) {
    return false;
  }
  return IsPrimary(ev, p);
}

// Particles of the given PDG code, either from the primary vertex (preFSI)
// or stable final state
inline std::vector<genie::GHepParticle *>
SelectParticles(genie::EventRecord const &ev, std::vector<int> const &pdgs,
                bool preFSI) {
  std::vector<genie::GHepParticle *> out;
  for (int i = 0; i < ev.GetEntries(); i++) {
    genie::GHepParticle *p = ev.Particle(i);
    bool pdg_match = false;
    for (int pdg : pdgs) {
      if (p->Pdg() == pdg) {
        pdg_match = true;
        break;
      }
    }
    if (!pdg_match) {
      continue;
    }
    if (preFSI ? IsPreFSIHadron(ev, p)
               : (p->Status() == genie::kIStStableFinalState)) {
      out.push_back(p);
    }
  }
  return out;
}

// Highest-momentum particle among the selection; nullptr if there is none
inline genie::GHepParticle *
LeadingParticle(std::vector<genie::GHepParticle *> const &particles) {
  genie::GHepParticle *leading = nullptr;
  for (genie::GHepParticle *p : particles) {
    if (!leading || p->P4()->Vect().Mag() > leading->P4()->Vect().Mag()) {
      leading = p;
    }
  }
  return leading;
}

// Adler angles; ported from NUISANCE FitUtils::CosThAdler and
// FitUtils::PhiAdler
inline double CosThAdler(TLorentzVector Pnu, TLorentzVector Pmu,
                         TLorentzVector Ppi, TLorentzVector Pprot) {
  TLorentzVector Pres = Pprot + Ppi;
  Pnu.Boost(-Pres.BoostVector());
  Pmu.Boost(-Pres.BoostVector());
  Ppi.Boost(-Pres.BoostVector());

  TVector3 zAxis =
      (Pnu.Vect() - Pmu.Vect()) * (1.0 / ((Pnu.Vect() - Pmu.Vect()).Mag()));

  return std::cos(Ppi.Vect().Angle(zAxis));
}

inline double PhiAdler(TLorentzVector Pnu, TLorentzVector Pmu,
                       TLorentzVector Ppi, TLorentzVector Pprot) {
  TLorentzVector Pres = Pprot + Ppi;
  Pnu.Boost(-Pres.BoostVector());
  Pmu.Boost(-Pres.BoostVector());
  Ppi.Boost(-Pres.BoostVector());

  TVector3 zAxis =
      (Pnu.Vect() - Pmu.Vect()) * (1.0 / ((Pnu.Vect() - Pmu.Vect()).Mag()));

  TVector3 yAxis = Pnu.Vect().Cross(Pmu.Vect());
  yAxis *= 1.0 / double(yAxis.Mag());

  TVector3 xAxis = yAxis.Cross(zAxis);
  xAxis *= 1.0 / double(xAxis.Mag());

  double x = Ppi.Vect().Dot(xAxis);
  double y = Ppi.Vect().Dot(yAxis);

  double newphi = std::atan2(y, x) * (180. / M_PI);
  if (newphi < 0.0) {
    newphi += 360.0;
  }
  return newphi;
}

// Adler angles need exactly one pion (any charge) and one nucleon, as in the
// NUISANCE flat tree. Returns {CosThAdler, PhiAdler}.
inline std::array<double, 2> AdlerAngles(genie::EventRecord const &ev,
                                         bool preFSI) {
  std::array<double, 2> out = {kUndefined, kUndefined};
  if (!HasFSLepton(ev)) {
    return out;
  }
  auto pions = SelectParticles(ev, {211, -211, 111}, preFSI);
  auto nucleons = SelectParticles(ev, {2212, 2112}, preFSI);
  if (pions.size() != 1 || nucleons.size() != 1) {
    return out;
  }
  TLorentzVector Pnu = ProbeP4(ev);
  TLorentzVector Plep = FSLeptonP4(ev);
  TLorentzVector Ppi = *pions[0]->P4();
  TLorentzVector PN = *nucleons[0]->P4();
  out[0] = CosThAdler(Pnu, Plep, Ppi, PN);
  out[1] = PhiAdler(Pnu, Plep, Ppi, PN);
  return out;
}

// PDG code of the final-state nucleon cluster of a GENIE MEC event
// (2000000200: nn, 2000000201: np, 2000000202: pp)
inline double MECFSClusterPDG(genie::EventRecord const &ev) {
  if (!ev.Summary()->ProcInfo().IsMEC()) {
    return kUndefined;
  }
  double pdg = kUndefined;
  for (int i = 0; i < ev.GetEntries(); i++) {
    int this_pdg = ev.Particle(i)->Pdg();
    if (this_pdg >= 2000000200 && this_pdg <= 2000000202) {
      pdg = this_pdg; // the last one is the final-state cluster
    }
  }
  return pdg;
}

// Pre-FSI 2p2h nucleon and lepton momenta in the reaction-plane frame, where
// the lepton transverse momentum lies along -y. Same convention as the
// SuSAv2 -> Valencia exclusive 2p2h BDT training (BDTReweight
// reaction_frame_2p2h_momenta):
//   {N1_px, N1_py, N1_pz, N2_px, N2_py, N2_pz, Lep_py, Lep_pz}
// N1 is the proton for an np pair and the higher-KE proton for a pp pair.
// Relies on the GENIE MEC record layout: 4 = lepton, 5 = final-state
// cluster, 6 and 7 = pre-FSI nucleons.
inline std::array<double, 8>
MECReactionPlaneMomenta(genie::EventRecord const &ev) {
  std::array<double, 8> out;
  out.fill(kUndefined);

  if (!ev.Summary()->ProcInfo().IsMEC() || ev.GetEntries() < 8) {
    return out;
  }

  genie::GHepParticle *lep = ev.Particle(4);
  genie::GHepParticle *cluster = ev.Particle(5);
  genie::GHepParticle *N1 = nullptr;
  genie::GHepParticle *N2 = nullptr;
  if (cluster->Pdg() == 2000000201) { // np
    bool first_is_p = ev.Particle(6)->Pdg() == 2212;
    N1 = first_is_p ? ev.Particle(6) : ev.Particle(7);
    N2 = first_is_p ? ev.Particle(7) : ev.Particle(6);
  } else if (cluster->Pdg() == 2000000202) { // pp
    bool first_is_leading = ev.Particle(6)->KinE() > ev.Particle(7)->KinE();
    N1 = first_is_leading ? ev.Particle(6) : ev.Particle(7);
    N2 = first_is_leading ? ev.Particle(7) : ev.Particle(6);
  } else {
    return out;
  }

  const double lep_px = lep->Px();
  const double lep_py = lep->Py();
  // new lepton py is negative of lepton transverse momentum magnitude
  const double lep_py_new = -std::sqrt(lep_px * lep_px + lep_py * lep_py);
  // unit vector y of reaction plane frame; x is y rotated clockwise by 90 deg
  const double y0 = lep_px / lep_py_new, y1 = lep_py / lep_py_new;
  const double x0 = -y1, x1 = y0;

  out[0] = N1->Px() * x0 + N1->Py() * x1;
  out[1] = N1->Px() * y0 + N1->Py() * y1;
  out[2] = N1->Pz();
  out[3] = N2->Px() * x0 + N2->Py() * x1;
  out[4] = N2->Px() * y0 + N2->Py() * y1;
  out[5] = N2->Pz();
  out[6] = lep_py_new;
  out[7] = lep->Pz();
  return out;
}

// ---------------------------------------------------------------------------
// Registry
// ---------------------------------------------------------------------------

inline std::map<std::string, EventVariable> BuildEventVariableRegistry() {
  std::map<std::string, EventVariable> reg;

  // Interaction
  reg["IsCC"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsWeakCC();
  };
  reg["IsNC"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsWeakNC();
  };
  reg["IsQE"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsQuasiElastic();
  };
  reg["IsMEC"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsMEC();
  };
  reg["IsRES"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsResonant();
  };
  reg["IsDIS"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsDeepInelastic();
  };
  reg["IsCOH"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ProcInfo().IsCoherentProduction();
  };
  reg["NeutMode"] = [](genie::EventRecord const &ev) -> double {
    return genie::utils::ghep::NeutReactionCode(&ev);
  };
  reg["ResonanceID"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->ExclTag().Resonance();
  };
  reg["ProbePDG"] = [](genie::EventRecord const &ev) -> double {
    return ev.Probe()->Pdg();
  };
  reg["TargetPDG"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->InitState().Tgt().Pdg();
  };
  reg["TargetA"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->InitState().Tgt().A();
  };
  reg["TargetZ"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->InitState().Tgt().Z();
  };
  // For MEC, this is the initial-state nucleon cluster PDG
  reg["HitNucleonPDG"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->InitState().Tgt().HitNucPdg();
  };
  reg["MECFSClusterPDG"] = MECFSClusterPDG;
  reg["HitNucleonP"] = [](genie::EventRecord const &ev) -> double {
    genie::GHepParticle *nucleon = ev.HitNucleon();
    return nucleon ? nucleon->P4()->Vect().Mag() : kUndefined;
  };
  reg["Emiss_GENIE"] = [](genie::EventRecord const &ev) -> double {
    genie::GHepParticle *nucleon = ev.HitNucleon();
    return nucleon ? nucleon->RemovalEnergy() : kUndefined;
  };

  // Neutrino and lepton
  reg["Enu"] = [](genie::EventRecord const &ev) -> double {
    return ProbeP4(ev).E();
  };
  reg["ELep"] = [](genie::EventRecord const &ev) -> double {
    return HasFSLepton(ev) ? FSLeptonP4(ev).E() : kUndefined;
  };
  reg["PLep"] = [](genie::EventRecord const &ev) -> double {
    return HasFSLepton(ev) ? FSLeptonP4(ev).Vect().Mag() : kUndefined;
  };
  reg["CosLep"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    return std::cos(ProbeP4(ev).Vect().Angle(FSLeptonP4(ev).Vect()));
  };
  // Lepton momentum along / perpendicular to the neutrino direction
  reg["PLepL"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    return FSLeptonP4(ev).Vect().Dot(ProbeP4(ev).Vect().Unit());
  };
  reg["PLepT"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    return FSLeptonP4(ev).Vect().Perp(ProbeP4(ev).Vect());
  };

  // Momentum transfer
  reg["q0"] = [](genie::EventRecord const &ev) -> double {
    return HasFSLepton(ev) ? MomentumTransfer(ev).E() : kUndefined;
  };
  reg["q3"] = [](genie::EventRecord const &ev) -> double {
    return HasFSLepton(ev) ? MomentumTransfer(ev).Vect().Mag() : kUndefined;
  };
  reg["Q2"] = [](genie::EventRecord const &ev) -> double {
    return HasFSLepton(ev) ? -MomentumTransfer(ev).Mag2() : kUndefined;
  };
  // W assuming the struck nucleon at rest
  reg["W_nuc_rest"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    TLorentzVector q = MomentumTransfer(ev);
    double W2 = kNucleonMass * kNucleonMass + 2. * kNucleonMass * q.E() +
                q.Mag2();
    return (W2 > 0.) ? std::sqrt(W2) : kUndefined;
  };
  reg["x_Bj"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    TLorentzVector q = MomentumTransfer(ev);
    return (q.E() > 0.) ? -q.Mag2() / (2. * kNucleonMass * q.E()) : kUndefined;
  };
  reg["y_inel"] = [](genie::EventRecord const &ev) -> double {
    if (!HasFSLepton(ev)) {
      return kUndefined;
    }
    return 1. - FSLeptonP4(ev).E() / ProbeP4(ev).E();
  };
  // Kinematics selected by GENIE
  reg["W_GENIE"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->Kine().W(true);
  };
  reg["Q2_GENIE"] = [](genie::EventRecord const &ev) -> double {
    return ev.Summary()->Kine().Q2(true);
  };

  // Missing energy / momentum, hadronic energy
  reg["Emiss"] = [](genie::EventRecord const &ev) -> double {
    return GetEmiss(ev, false);
  };
  reg["Emiss_preFSI"] = [](genie::EventRecord const &ev) -> double {
    return GetEmiss(ev, true);
  };
  reg["pmiss"] = [](genie::EventRecord const &ev) -> double {
    return GetPmiss(ev, false);
  };
  reg["pmiss_preFSI"] = [](genie::EventRecord const &ev) -> double {
    return GetPmiss(ev, true);
  };
  reg["EAvail"] = [](genie::EventRecord const &ev) -> double {
    return GetErecoil_MINERvA_LowRecoil(ev);
  };

  // Hadron multiplicity and leading-hadron kinematics, post-FSI and pre-FSI:
  //   N<X>: multiplicity, P<X>: |p|, T<X>: kinetic energy, E<X>: energy,
  //   Cos<X>: cosine of the angle to the neutrino
  // of the highest-momentum <X>, where <X> is one of P, N, PiP, PiM, Pi0.
  // E.g. "NPiP", "preFSI_PP", "preFSI_CosP".
  std::vector<std::pair<std::string, int>> const species = {
      {"P", 2212}, {"N", 2112}, {"PiP", 211}, {"PiM", -211}, {"Pi0", 111}};
  for (bool preFSI : {false, true}) {
    std::string const prefix = preFSI ? "preFSI_" : "";
    for (auto const &[name, pdg] : species) {
      reg[prefix + "N" + name] = [pdg = pdg, preFSI](
                                     genie::EventRecord const &ev) -> double {
        return SelectParticles(ev, {pdg}, preFSI).size();
      };
      reg[prefix + "P" + name] = [pdg = pdg, preFSI](
                                     genie::EventRecord const &ev) -> double {
        auto p = LeadingParticle(SelectParticles(ev, {pdg}, preFSI));
        return p ? p->P4()->Vect().Mag() : kUndefined;
      };
      reg[prefix + "T" + name] = [pdg = pdg, preFSI](
                                     genie::EventRecord const &ev) -> double {
        auto p = LeadingParticle(SelectParticles(ev, {pdg}, preFSI));
        return p ? p->KinE() : kUndefined;
      };
      reg[prefix + "E" + name] = [pdg = pdg, preFSI](
                                     genie::EventRecord const &ev) -> double {
        auto p = LeadingParticle(SelectParticles(ev, {pdg}, preFSI));
        return p ? p->E() : kUndefined;
      };
      reg[prefix + "Cos" + name] = [pdg = pdg, preFSI](
                                       genie::EventRecord const &ev) -> double {
        auto p = LeadingParticle(SelectParticles(ev, {pdg}, preFSI));
        return p ? std::cos(p->P4()->Vect().Angle(ProbeP4(ev).Vect()))
                 : kUndefined;
      };
    }

    // Adler angles of the single-pion, single-nucleon system
    reg[prefix + "CosThAdler"] = [preFSI](
                                     genie::EventRecord const &ev) -> double {
      return AdlerAngles(ev, preFSI)[0];
    };
    reg[prefix + "PhiAdler"] = [preFSI](
                                   genie::EventRecord const &ev) -> double {
      return AdlerAngles(ev, preFSI)[1];
    };
  }

  // 2p2h reaction-plane-frame momenta
  std::vector<std::string> const rp_names = {
      "MEC_RP_N1_px", "MEC_RP_N1_py", "MEC_RP_N1_pz",  "MEC_RP_N2_px",
      "MEC_RP_N2_py", "MEC_RP_N2_pz", "MEC_RP_Lep_py", "MEC_RP_Lep_pz"};
  for (size_t i = 0; i < rp_names.size(); i++) {
    reg[rp_names[i]] = [i](genie::EventRecord const &ev) -> double {
      return MECReactionPlaneMomenta(ev)[i];
    };
  }

  return reg;
}

inline std::map<std::string, EventVariable> const &GetEventVariableRegistry() {
  static std::map<std::string, EventVariable> const reg =
      BuildEventVariableRegistry();
  return reg;
}

inline std::vector<std::string> GetAvailableEventVariables() {
  std::vector<std::string> names;
  for (auto const &[name, fn] : GetEventVariableRegistry()) {
    names.push_back(name);
  }
  return names;
}

inline EventVariable GetEventVariable(std::string const &name) {
  auto const &reg = GetEventVariableRegistry();
  auto it = reg.find(name);
  if (it == reg.end()) {
    std::stringstream ss;
    for (auto const &avail : GetAvailableEventVariables()) {
      ss << " " << avail;
    }
    throw unknown_event_variable()
        << "[ERROR]: Unknown event variable " << std::quoted(name)
        << ". Available variables:" << ss.str();
  }
  return it->second;
}

// ---------------------------------------------------------------------------
// Cuts
// ---------------------------------------------------------------------------

// A cut "<variable> <op> <value>", e.g. "Enu < 10" or "NeutMode == 11",
// where <op> is one of ==, !=, <, <=, >, >=
struct EventCut {
  std::string expression;
  EventVariable variable;
  std::string op;
  double value;

  bool Pass(genie::EventRecord const &ev) const {
    double x = variable(ev);
    if (op == "==") return x == value;
    if (op == "!=") return x != value;
    if (op == "<") return x < value;
    if (op == "<=") return x <= value;
    if (op == ">") return x > value;
    return x >= value; // ">="
  }
};

inline EventCut ParseEventCut(std::string const &expression) {
  size_t op_pos = expression.find_first_of("=!<>");
  if (op_pos == std::string::npos) {
    throw invalid_event_cut()
        << "[ERROR]: No comparison operator in cut "
        << std::quoted(expression);
  }
  size_t op_len = (op_pos + 1 < expression.size() &&
                   expression[op_pos + 1] == '=')
                      ? 2
                      : 1;
  std::string op = expression.substr(op_pos, op_len);
  if (op == "=" || op == "!") {
    throw invalid_event_cut() << "[ERROR]: Invalid operator " << std::quoted(op)
                              << " in cut " << std::quoted(expression)
                              << "; use one of ==, !=, <, <=, >, >=";
  }

  std::string lhs = expression.substr(0, op_pos);
  std::string rhs = expression.substr(op_pos + op_len);
  systtools::trim(lhs);
  systtools::trim(rhs);

  double value = 0.;
  size_t n_parsed = 0;
  try {
    value = std::stod(rhs, &n_parsed);
  } catch (std::exception const &) {
    n_parsed = 0;
  }
  if (rhs.empty() || n_parsed != rhs.size()) {
    throw invalid_event_cut() << "[ERROR]: Could not parse value "
                              << std::quoted(rhs) << " in cut "
                              << std::quoted(expression);
  }

  return EventCut{expression, GetEventVariable(lhs), op, value};
}

} // namespace evvar
} // namespace nusyst
