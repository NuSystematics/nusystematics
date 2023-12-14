#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/MINERvARPAq0q3_ReWeight.hh"
#include "nusystematics/responsecalculators/MINERvAq0q3Weighting_data.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

///
///
/// Can apply the MINERvA RPA tune and two discrete 'tweaks', the 2p2h gaussian
/// enhancement tune and 3 discrete 'tweaks', can also apply the 2p2h gaussian
/// with arbitrary parameters.
class MINERvAq0q3Weighting : public nusyst::IGENIESystProvider_tool {

  enum class param_t {
    kGaussNorm = 0,
    kGaussMeanQ0,
    kGaussMeanQ3,
    kGaussSigmaQ0,
    kGaussSigmaQ3,
    kGaussCorrelation,
    kGaussResponse,
    kMINERvARPA,
    kMINERvA2p2h,
    kMINERvA2p2h_CV,
    kMINERvA2p2h_NN,
    kMINERvA2p2h_np,
    kMINERvA2p2h_QE,
  };

  struct param_t_name {
    // The parameter prettyName expected to be found in the input metadata
    std::string name;
    // The local enum id of the parameter
    param_t lid;
  };

  std::unique_ptr<nusyst::MINERvARPAq0q3_ReWeight> RPATemplateReweighter;
  std::map<param_t, size_t> ConfiguredParameters;

public:
  explicit MINERvAq0q3Weighting(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  double GetMINERvARPATuneWeight(double val, double q0, double q3);
  double GetMINERvA2p2hTuneEnhancement(int val, double q0, double q3,
                                       nusyst::QELikeTarget_t QELTarget);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~MINERvAq0q3Weighting();

private:
  fhicl::ParameterSet tool_options;
  std::pair<double, double> MEC_LimitWeights;

  void InitValidTree();

  std::vector<double> vals_2p2hTotal, vals_2p2hCV, vals_2p2hNN, vals_2p2hnp,
      vals_2p2hQE;

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  bool parameter_per_2p2h_universe;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget, nRPA_weights, nMEC_weights;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
  std::vector<double> RPA_weights;
  std::vector<double> MEC_weights;
};
