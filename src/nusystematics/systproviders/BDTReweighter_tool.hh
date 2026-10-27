// File: BDTReweighter_tool.hh
//
// Generic reweighter using BDT (hep_ml GBReweighter) models exported to json
// with the BDTReweight python module.
// The json model file, the ordered list of input variables, and the event
// selection are all given in the tool configuration; variable names are
// resolved with nusyst::evvar (utility/GENIEEventVariables.hh).

#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/utility/BDTReweighter_json.hh"
#include "nusystematics/utility/GENIEEventVariables.hh"

#include <memory>
#include <string>
#include <vector>

class BDTReweighter : public nusyst::IGENIESystProvider_tool {

public:
  explicit BDTReweighter(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  // Extra options are saved in tool_options
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString() { return "BDTReweighter"; };

  ~BDTReweighter();

private:
  // One BDT model applied to the events passing its selection
  struct BDTModel {
    std::string name;
    std::unique_ptr<BDTReweight::JSONReweighter> reweighter;
    std::vector<std::string> variable_names;
    std::vector<nusyst::evvar::EventVariable> variables;
    std::vector<nusyst::evvar::EventCut> selection;
    double scale_factor;
    double weight_min;
    double weight_max;

    bool PassSelection(genie::EventRecord const &ev) const;
    double GetWeight(genie::EventRecord const &ev) const;
  };

  fhicl::ParameterSet tool_options;
  std::string dial_name;
  size_t pidx_Dial;
  std::vector<BDTModel> models;
};
