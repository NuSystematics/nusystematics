// File: BDTReweighter_tool.cc

#include "nusystematics/systproviders/BDTReweighter_tool.hh"
#include "nusystematics/utility/exceptions.hh"
#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace systtools;
using namespace nusyst;

BDTReweighter::BDTReweighter(YAML::Node const &params)
    : IGENIESystProvider_tool(params),
      pidx_Dial(systtools::kParamUnhandled<size_t>) {}

SystMetaData BDTReweighter::BuildSystMetaData(YAML::Node const &cfg,
                                              paramId_t firstId) {

  std::cout << "[BDTReweighter::BuildSystMetaData] Called" << std::endl;
  SystMetaData smd;

  // Name of the dial; must be unique when this tool is used more than once
  dial_name = cfg["dial_name"] ? cfg["dial_name"].as<std::string>()
                               : "BDTReweight";
  tool_options["dial_name"] = dial_name;

  systtools::SystParamHeader phdr;
  if (ParseYAMLToolConfigurationParameter(cfg, dial_name, phdr, firstId)) {
    printf("[BDTReweighter::BuildSystMetaData] %s is found from ToolConfig\n",
           dial_name.c_str());
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  } else {
    std::cout << "[BDTReweighter::BuildSystMetaData] No dial is set"
              << std::endl;
  }

  if (!cfg["reweighters"] || !cfg["reweighters"].IsSequence()) {
    throw std::runtime_error(
        "BDTReweighter: reweighters must be a list of BDT model configs");
  }
  tool_options["reweighters"] = cfg["reweighters"];

  return smd;
}

bool BDTReweighter::SetupResponseCalculator(YAML::Node const &tool_options) {

  std::cout << "[BDTReweighter::SetupResponseCalculator] Called" << std::endl;
  systtools::SystMetaData const &md = GetSystMetaData();

  dial_name = tool_options["dial_name"]
                  ? tool_options["dial_name"].as<std::string>()
                  : "BDTReweight";
  if (HasParam(md, dial_name)) {
    pidx_Dial = GetParamIndex(md, dial_name);
  }

  models.clear();
  for (YAML::Node const &model_cfg : tool_options["reweighters"]) {
    BDTModel model;
    model.name = model_cfg["name"] ? model_cfg["name"].as<std::string>()
                                   : "BDT_" + std::to_string(models.size());
    std::string const prefix =
        "[BDTReweighter::SetupResponseCalculator][" + model.name + "] ";

    std::string json_path = model_cfg["json_path"]
                                ? model_cfg["json_path"].as<std::string>()
                                : "";
    if (json_path.empty()) {
      throw std::runtime_error("BDTReweighter: json_path is empty for " +
                               model.name);
    }
    std::cout << prefix << "json_path: " << json_path << std::endl;
    model.reweighter = std::make_unique<BDTReweight::JSONReweighter>(json_path);

    // Input variables, in the same order as used in the training
    if (!model_cfg["variables"]) {
      throw std::runtime_error("BDTReweighter: variables is not set for " +
                               model.name);
    }
    model.variable_names =
        model_cfg["variables"].as<std::vector<std::string>>();
    if (model.variable_names.size() != model.reweighter->NFeatures()) {
      throw std::runtime_error(
          "BDTReweighter: " + model.name + " has " +
          std::to_string(model.variable_names.size()) +
          " variables, but the BDT model expects " +
          std::to_string(model.reweighter->NFeatures()));
    }
    std::cout << prefix << "variables:";
    for (std::string const &vname : model.variable_names) {
      model.variables.push_back(evvar::GetEventVariable(vname));
      std::cout << " " << vname;
    }
    std::cout << std::endl;

    // Event selection; all cuts must pass. Events not selected by any model
    // get unit weight.
    std::vector<std::string> cuts =
        model_cfg["selection"]
            ? model_cfg["selection"].as<std::vector<std::string>>()
            : std::vector<std::string>{};
    for (std::string const &cut : cuts) {
      model.selection.push_back(evvar::ParseEventCut(cut));
      std::cout << prefix << "selection: " << cut << std::endl;
    }

    // BDT reweight is shape-only; the normalization is set by scale_factor
    model.scale_factor = model_cfg["scale_factor"]
                             ? model_cfg["scale_factor"].as<double>()
                             : 1.0;
    // Optional clipping of the BDT weight, applied before scale_factor
    model.weight_min = model_cfg["weight_min"]
                           ? model_cfg["weight_min"].as<double>()
                           : 0.;
    model.weight_max = model_cfg["weight_max"]
                           ? model_cfg["weight_max"].as<double>()
                           : std::numeric_limits<double>::infinity();
    std::cout << prefix << "scale_factor: " << model.scale_factor
              << ", weight_min: " << model.weight_min
              << ", weight_max: " << model.weight_max << std::endl;

    models.push_back(std::move(model));
  }

  return true;
}

bool BDTReweighter::BDTModel::PassSelection(
    genie::EventRecord const &ev) const {
  for (evvar::EventCut const &cut : selection) {
    if (!cut.Pass(ev)) {
      return false;
    }
  }
  return true;
}

double
BDTReweighter::BDTModel::GetWeight(genie::EventRecord const &ev) const {
  std::vector<double> features;
  features.reserve(variables.size());
  for (evvar::EventVariable const &var : variables) {
    features.push_back(var(ev));
  }
  double weight = reweighter->PredictWeight(features);
  weight = std::clamp(weight, weight_min, weight_max);
  return weight * scale_factor;
}

event_unit_response_t
BDTReweighter::GetEventResponse(genie::EventRecord const &ev) {

  if (pidx_Dial == systtools::kParamUnhandled<size_t>) {
    return this->GetDefaultEventResponse();
  }

  // The first model whose selection passes is used
  BDTModel const *model = nullptr;
  for (BDTModel const &m : models) {
    if (m.PassSelection(ev)) {
      model = &m;
      break;
    }
  }
  if (!model) {
    return this->GetDefaultEventResponse();
  }

  double weight = model->GetWeight(ev);

  systtools::event_unit_response_t resp;
  systtools::SystMetaData const &md = GetSystMetaData();
  resp.push_back({md[pidx_Dial].systParamId, {}});
  for (double var : md[pidx_Dial].paramVariations) {
    // var = 0: CV (unit weight), var = 1: fully reweighted
    resp.back().responses.push_back(1 + var * (weight - 1));
  }

  return resp;
}

BDTReweighter::~BDTReweighter() {}
