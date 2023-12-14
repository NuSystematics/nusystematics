#pragma once

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

// #define TemplateResponseCalculatorBase_DEBUG

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(no_responses_loaded);
NEW_SYSTTOOLS_EXCEPT(incompatible_number_of_bins);
NEW_SYSTTOOLS_EXCEPT(bad_value_ordering);

template <size_t NDims, bool Continuous = true, size_t PolyResponseOrder = 5>
class TemplateResponseCalculatorBase {

protected:
  std::vector<systtools::PolyResponse<PolyResponseOrder>>
      InterpolatedBinResponses;
  std::map<double, std::unique_ptr<typename THType<NDims>::type>>
      BinnedResponses;

  void ValidateInputHistograms();
  void BuildInterpolatedResponses();

public:
  static size_t const NDimensions = NDims;
  TemplateResponseCalculatorBase();
  TemplateResponseCalculatorBase(TemplateResponseCalculatorBase &&other)
      : InterpolatedBinResponses(std::move(other.InterpolatedBinResponses)),
        BinnedResponses(std::move(other.BinnedResponses)) {}

  /// Reads and loads input fhicl
  ///
  /// Expected fhicl like:
  ///  use_FW_SEARCH_PATH: true # If enabled, will search for files in
  ///  in the search path. Only read for ART jobs.
  ///  input_file: "file.root" # Optional default root file name for this
  ///                          # parameter's inputs
  ///    inputs: [
  ///      { value: 0
  ///        input_file: "file.root" # Optional if the less-specific is
  ///                                # supplied
  ///        input_hist: "histo_name"
  ///      }
  ///    ]
  ///  }
  void LoadInputHistograms(fhicl::ParameterSet const &ps);

  typedef Int_t bin_it_t;

  virtual bin_it_t GetBin(std::array<double, NDims> const &) const;

  virtual std::string GetCalculatorName() const = 0;
  
  // virtual destructor required
  virtual ~TemplateResponseCalculatorBase() = default;

  /// Extra template parameters required for SINFAE
  template <bool IsCont = Continuous>
  typename std::enable_if<IsCont, double>::type
  GetVariation(double, typename std::enable_if<IsCont, bin_it_t>::type) const;
  /// Extra template parameters required for SINFAE
  template <bool IsCont = Continuous>
  typename std::enable_if<!IsCont, double>::type
  GetVariation(double, typename std::enable_if<!IsCont, bin_it_t>::type) const;

  double GetVariation(double val,
                      std::array<double, NDims> const &kinematics) const;

  std::vector<double> GetValidVariations() const;
  bool IsValidVariation(double val) const;
};

//*********************** Implementations

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
void TemplateResponseCalculatorBase<
    NDims, Continuous, PolyResponseOrder>::ValidateInputHistograms() {
  if (!BinnedResponses.size()) {
    throw no_responses_loaded()
        << "[ERROR]: Expected to find some responses, but found none.";
  }
  if (Continuous && (BinnedResponses.size() == 1)) {
    throw no_responses_loaded()
        << "[ERROR]: For continuous template parameter, only a single set of "
           "responses was loaded, require at least two parameter values for "
           "continuous response.";
  }
  size_t NBins = THType<NDims>::GetNbins(BinnedResponses.begin()->second);
  for (auto &val_resp : BinnedResponses) {
    if (THType<NDims>::GetNbins(val_resp.second) != NBins) {
      throw incompatible_number_of_bins()
          << "[ERROR]: The first histogram at parameter value "
          << BinnedResponses.begin()->first << " has a response in " << NBins
          << " bins, at parameter value " << val_resp.first << " found "
          << THType<NDims>::GetNbins(val_resp.second) << " bins.";
    }
  }
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
void TemplateResponseCalculatorBase<NDims, Continuous, PolyResponseOrder>::
    LoadInputHistograms(fhicl::ParameterSet const &ps) {

  std::string const &default_root_file = ps.get<std::string>("input_file", "");

  for (fhicl::ParameterSet const &val_config :
       ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
    double pval = val_config.get<double>("value");
    std::string input_file =
        val_config.get<std::string>("input_file", default_root_file);
    std::string input_hist = val_config.get<std::string>("input_hist");

    BinnedResponses[pval] = std::unique_ptr<typename THType<NDims>::type>(
        GetHistogram<typename THType<NDims>::type>(input_file, input_hist));
  }

  ValidateInputHistograms();
  if (Continuous) {
    BuildInterpolatedResponses();
  }
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
typename TemplateResponseCalculatorBase<NDims, Continuous,
                                        PolyResponseOrder>::bin_it_t
TemplateResponseCalculatorBase<NDims, Continuous, PolyResponseOrder>::GetBin(
    std::array<double, NDims> const &vals) const {
  return THType<NDims>::GetBin(BinnedResponses.begin()->second.get(), vals);
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
void TemplateResponseCalculatorBase<
    NDims, Continuous, PolyResponseOrder>::BuildInterpolatedResponses() {
  std::vector<double> xvals;
  std::vector<double> yvals_dummy;
  std::vector<double> yvals;
  for (auto const &var : BinnedResponses) {
    if (var.first < xvals.back()) {
      throw bad_value_ordering()
          << "[ERROR]: When precalculating response functions, found value "
             "specification for "
          << var.first << ", but the previous value was " << xvals.back()
          << "; as these are stored in a std::map, this is problematic...";
    }
    xvals.push_back(var.first);
    yvals_dummy.push_back(1);
  }

  size_t NBins = THType<NDims>::GetNbins(BinnedResponses.begin()->second, true);
  for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
    yvals.clear();
    for (auto const &var : BinnedResponses) {
      if (THType<NDims>::IsFlowBin(var.second, bi_it)) {
        yvals = yvals_dummy;
        break;
      }
      yvals.push_back(var.second->GetBinContent(bi_it));
    }
    InterpolatedBinResponses.emplace_back(xvals, yvals);
  }
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
TemplateResponseCalculatorBase<
    NDims, Continuous, PolyResponseOrder>::TemplateResponseCalculatorBase() {}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
template <bool IsCont>
typename std::enable_if<IsCont, double>::type
TemplateResponseCalculatorBase<NDims, Continuous, PolyResponseOrder>::
    GetVariation(double val,
                 typename std::enable_if<IsCont, bin_it_t>::type bin) const {
  return InterpolatedBinResponses[bin]->eval(val);
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
template <bool IsCont>
typename std::enable_if<!IsCont, double>::type
TemplateResponseCalculatorBase<NDims, Continuous, PolyResponseOrder>::
    GetVariation(double val,
                 typename std::enable_if<!IsCont, bin_it_t>::type bin) const {

  if (bin == kBinOutsideRange) {
    return 1;
  }

  for (auto const &resp : BinnedResponses) {
    if (fabs(val - resp.first) <
        (std::numeric_limits<double>::epsilon() * 1E4)) {
#ifdef TemplateResponseCalculatorBase_DEBUG
      std::cout << "[INFO]: Getting bin content for bin: " << bin
                << " at value: " << val << " = "
                << resp.second->GetBinContent(bin)
                << " from hist: " << resp.second->GetName() << std::endl;
#endif
      return resp.second->GetBinContent(bin);
    }
  }

  std::stringstream ss("");
  ss << "[";
  for (auto const &v : GetValidVariations()) {
    ss << v << ", ";
  }
  std::string valid_vals = ss.str();
  valid_vals = valid_vals.substr(0, valid_vals.size() - 2) + " ]";
  throw systtools::invalid_parameter_value()
      << "[ERROR]: Invalid parameter value, " << val
      << " used for template response " << GetCalculatorName()
      << ", configured values: " << valid_vals;
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
double TemplateResponseCalculatorBase<NDims, Continuous, PolyResponseOrder>::
    GetVariation(double val,
                 std::array<double, NDims> const &kinematics) const {
  return GetVariation(val, GetBin(kinematics));
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
std::vector<double>
TemplateResponseCalculatorBase<NDims, Continuous,
                               PolyResponseOrder>::GetValidVariations() const {
  std::vector<double> vals;
  if (Continuous) {
    std::pair<double, double> minmax{std::numeric_limits<double>::max(),
                                     -std::numeric_limits<double>::max()};
    for (auto const &var : BinnedResponses) {
      minmax.first = std::min(minmax.first, var.first);
      minmax.second = std::max(minmax.second, var.first);
    }
    vals.push_back(minmax.first);
    vals.push_back(minmax.second);
  } else {
    for (auto const &var : BinnedResponses) {
      vals.push_back(var.first);
    }
  }
  return vals;
}

template <size_t NDims, bool Continuous, size_t PolyResponseOrder>
bool TemplateResponseCalculatorBase<
    NDims, Continuous, PolyResponseOrder>::IsValidVariation(double val) const {
  std::vector<double> const &valvar = GetValidVariations();
  if (Continuous) {
    return (val > valvar[0]) && (val < valvar[1]);
  } else {
    for (double v : valvar) {
      if (fabs(v - val) < (std::numeric_limits<double>::epsilon() * 1E4)) {
        return true;
      }
    }
    return false;
  }
}

typedef TemplateResponseCalculatorBase<2, false> TemplateResponse2DDiscrete;
typedef TemplateResponseCalculatorBase<3, false> TemplateResponse3DDiscrete;

} // namespace nusyst
