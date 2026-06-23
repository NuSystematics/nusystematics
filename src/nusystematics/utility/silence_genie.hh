#pragma once

// Helpers to make the nusystematics CLI tools' output user-readable by
// suppressing the various sources of GENIE / ROOT / provider chatter that
// otherwise dominate the screen:
//
//   * GENIE's mandatory ASCII banner — printed unconditionally via cout by
//     genie::Messenger::Instance() the first time it's invoked. Not gated
//     by any log priority or env var; the only way to suppress is to swap
//     cout's stream buffer while the first call happens.
//
//   * GENIE log4cpp INFO/NOTICE messages (RunOpt::BuildTune,
//     NievesQELCCPXSec::LoadConfig, TransverseEnhancementFFModel, ...).
//     Suppressed by Messenger_whisper.xml, which raises priorities to WARN+.
//
//   * nusyst provider [INFO]: / [LOUD]: lines emitted during construction
//     (ZExpPCAWeighter, IGENIESystProvider, ...). Direct cout writes; only
//     a stream-buffer swap silences them.
//
//   * ROOT INFO / WARN messages routed through gErrorIgnoreLevel — the
//     "Replacing existing TGraph2D" warning is the canonical example.
//
// Typical usage at the top of a CLI tool's main():
//
//   nusyst::quiet::SetGlobalQuiet();           // ROOT + GENIE log priorities
//   {
//     nusyst::quiet::StdoutSink _quiet;        // hide banner + cout chatter
//     genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
//         "Messenger_whisper.xml");
//     response_helper phh(cliopts::fclname);   // noisy GENIE/provider init
//   }
//
// Once the StdoutSink scope ends, cout is restored and normal tool output
// (progress lines, tables, etc.) prints as usual. cerr is never touched, so
// genuine errors and the [OK]/[FAIL]/[SKIP] status lines (which the config
// tool routes to cerr) keep flowing.

#include "TError.h"

#include "Framework/Messenger/Messenger.h"

#include <fstream>
#include <iostream>

namespace nusyst::quiet {

// RAII helper: while live, std::cout is redirected to /dev/null. Non-copyable.
class StdoutSink {
 public:
  StdoutSink() : devnull_("/dev/null"), orig_(nullptr) {
    std::cout.flush();
    orig_ = std::cout.rdbuf(devnull_.rdbuf());
  }
  ~StdoutSink() {
    std::cout.flush();
    if (orig_) std::cout.rdbuf(orig_);
  }
  StdoutSink(const StdoutSink &) = delete;
  StdoutSink &operator=(const StdoutSink &) = delete;
 private:
  std::ofstream devnull_;
  std::streambuf *orig_;
};

// One-shot global mute for ROOT diagnostics. Idempotent.
// (Does NOT touch GENIE's Messenger — its Instance() call triggers the
// banner, so callers should construct a StdoutSink around the first GENIE
// touch and explicitly call SetPrioritiesFromXmlFile inside that scope.)
inline void SetGlobalQuiet() {
  static bool done = false;
  if (done) return;
  gErrorIgnoreLevel = kError;
  done = true;
}

}  // namespace nusyst::quiet
