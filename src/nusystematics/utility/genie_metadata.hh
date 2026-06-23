#pragma once

// Inspect a GENIE GHEP file's metadata to extract the generator version and
// tune that produced it. Used to surface GENIE-version-skew issues at runtime
// rather than failing mysteriously inside a Jacobian transform (see the long
// debugging story of `kPSTAfE`-on-QE events from GENIE 3.04 samples).
//
// A populated GHEP file carries a `genie::NtpMCTreeHeader` TKey with two
// TObjStrings: `cvstag` (intended to record the generator's version string)
// and `tune`. The cvstag is, in practice, often left at the placeholder
// "999.999.999" when GENIE wasn't built with the version-stamping macros
// enabled — we treat that as "could not detect".

#include "Framework/Ntuple/NtpMCTreeHeader.h"

#include "TFile.h"

#include <iostream>
#include <optional>
#include <string>

namespace nusyst::metadata {

struct GENIESampleInfo {
  std::string cvstag; // raw value from the file ("999.999.999" → placeholder)
  std::string tune;
  bool cvstag_is_real() const {
    if (cvstag.empty()) return false;
    if (cvstag == "999.999.999") return false; // GENIE default placeholder
    return true;
  }
};

inline std::optional<GENIESampleInfo> ReadSampleInfo(const std::string &path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) return std::nullopt;
  auto *hdr = dynamic_cast<genie::NtpMCTreeHeader *>(f->Get("header"));
  if (!hdr) return std::nullopt;
  GENIESampleInfo info;
  info.cvstag = hdr->cvstag.GetString().Data();
  info.tune = hdr->tune.GetString().Data();
  return info;
}

// Convenience: print a one-line report to `os`. The caller hands in the input
// path that's already been resolved (i.e. a single .root file, not a TChain
// wildcard descriptor — passing a wildcard descriptor will simply fail to
// open and we'll fall through to "could not detect").
inline void ReportSampleInfo(const std::string &path,
                              std::ostream &os = std::cout) {
  auto info = ReadSampleInfo(path);
  if (!info) {
    os << "[INFO]: GENIE sample metadata: could not detect "
       << "(no NtpMCTreeHeader in " << path << ")" << std::endl;
    return;
  }
  if (info->cvstag_is_real()) {
    os << "[INFO]: GENIE sample metadata: version=" << info->cvstag
       << ", tune=" << info->tune << std::endl;
  } else {
    os << "[INFO]: GENIE sample metadata: version could not be detected "
       << "(cvstag=\"" << info->cvstag << "\")"
       << ", tune=" << info->tune << std::endl;
  }
}

} // namespace nusyst::metadata
