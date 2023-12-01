#pragma once

#include "systematicstools/utility/ROOTUtility.hh"

#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

#include "systematicstools/utility/string_parsers.hh"

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(non_contiguous_enu_range);

/// Similar to TemplateResponseCaluclatorBase but handles a Enu 'axis' manually
/// to allow non-uniform binning in the other dimensions.
///
/// Must be specialized with a TemplateResponseCalculator subclass that knows
/// how to search for bins, see MKSinglePiTemplate_ReWeight.hh for an example
template <class TRC> class EnuBinnedTemplateResponseCalculator {
public:
  typedef Int_t enu_bin_it_t;

private:
  static std::string StringifyNumberToOneDP(double number) {
    std::stringstream ss("");
    if (fabs(number - double(int(number))) < 1E-8) {
      ss << int(number);
    } else {
      ss << std::fixed << std::setprecision(2) << number;
    }
    return ss.str();
  }

  std::vector<double> EnuBinning;
  std::vector<TRC> EnuResponses;

  /// Reads and loads input fhicl
  ///
  /// Expected fhicl like:
  /// <inputs_key>: {
  ///    input_file_pattern: "example.%EGeV_%V.root" # optional, %E replaced
  ///                                                # with enu bin center in
  ///                                                # GeV to one decimal
  ///                                                # place, %V replace with
  ///                                                # parameter value to one
  ///                                                # decimal place.
  ///    input_hist_pattern: "hist_%EGeV_%V" # optional, %E replaced with enu
  ///                                        # bin center in GeV to one decimal
  ///                                        # place, %V replace with parameter
  ///                                        # value to one  decimal place.
  ///    e_uniform: [<nbins>,<e_low>,<e_high>] # optional
  //     param_values: [<val1>,<val2>,...] # optional
  ///    e_stops: [ {
  ///      enu_range: [1.5,2.5] # optional
  ///      value_uniform: [<nvalues>,<v_min>,<v_max>] # optional
  ///      input_file_pattern: example.1GeV_%V.root # optional
  ///      input_hist_pattern: myhist_2GeV_%V # optional
  ///      param_values: [ {
  ///          value: 1 # optional
  ///          input_file: example.1GeV_2.root # optional
  ///          input_hist: myhist_2GeV_2 # optional
  ///        }
  ///      ] # optional if <inputs_key>.param_values has been specified
  ///      }
  ///    ] # optional if all of input_file_pattern, input_hist_pattern,
  ///      # e_uniform, and param_values are specified
  /// }
  void LoadInputHistograms(fhicl::ParameterSet const &ps) {
    bool uniform_enu = false;
    bool consistent_param_values = false;

    std::vector<double> param_values;

    if (ps.has_key("e_uniform")) {
      std::tuple<size_t, double, double> e_uniform_descriptor =
          ps.get<std::tuple<size_t, double, double>>("e_uniform");

      double step = (std::get<2>(e_uniform_descriptor) -
                     std::get<1>(e_uniform_descriptor)) /
                    double(std::get<0>(e_uniform_descriptor));

      for (size_t i = 0; i <= std::get<0>(e_uniform_descriptor); ++i) {
        EnuBinning.push_back(std::get<1>(e_uniform_descriptor) + i * step);
      }
      uniform_enu = true;
    }

    if (ps.has_key("param_values")) {
      consistent_param_values = true;
      param_values = ps.get<std::vector<double>>("param_values");
    }

    std::string input_file_pattern =
        ps.get<std::string>("input_file_pattern", "");
    std::string input_hist_pattern =
        ps.get<std::string>("input_hist_pattern", "");

    if (!ps.has_key("e_stops")) { // all specified by patterns
      if (!consistent_param_values || !uniform_enu) {
        throw;
      }
      std::vector<fhicl::ParameterSet> value_descriptors;
      for (size_t e_stop_ctr = 0; e_stop_ctr < (EnuBinning.size() - 1);
           ++e_stop_ctr) {
        std::pair<double, double> enu_range =
            std::make_pair(EnuBinning[e_stop_ctr], EnuBinning[e_stop_ctr + 1]);
        double enu_mid = (enu_range.second + enu_range.first) / 2.0;

        std::string input_file_pattern_stop = systtools::str_replace(
            input_file_pattern, "%E", StringifyNumberToOneDP(enu_mid));
        std::string input_hist_pattern_stop = systtools::str_replace(
            input_hist_pattern, "%E", StringifyNumberToOneDP(enu_mid));

        for (size_t pval_ctr = 0; pval_ctr < param_values.size(); ++pval_ctr) {
          std::string input_file = systtools::str_replace(
              input_file_pattern_stop, "%V",
              StringifyNumberToOneDP(param_values[pval_ctr]));
          std::string input_hist = systtools::str_replace(
              input_hist_pattern_stop, "%V",
              StringifyNumberToOneDP(param_values[pval_ctr]));

          fhicl::ParameterSet value_descriptor;
          value_descriptor.put("value", param_values[pval_ctr]);
          value_descriptor.put("input_file", input_file);
          value_descriptor.put("input_hist", input_hist);
          value_descriptors.push_back(std::move(value_descriptor));
        }
        fhicl::ParameterSet estop_descriptor;
        estop_descriptor.put("inputs", value_descriptors);
        EnuResponses.emplace_back();
        EnuResponses.back().LoadInputHistograms(estop_descriptor);
      }
      return;
    }

    size_t e_stop_ctr = 0;
    for (fhicl::ParameterSet const &e_stop :
         ps.get<std::vector<fhicl::ParameterSet>>("e_stops")) {

      fhicl::ParameterSet estop_descriptor;
      std::pair<double, double> enu_range;

      if (!e_stop.has_key("enu_range") && uniform_enu) {
        if (EnuBinning.size() >= (e_stop_ctr + 1)) {
          throw incompatible_number_of_bins()
              << "[ERROR]: Evaluating e_stop #" << e_stop_ctr
              << " with uniform enu range, but only have " << EnuBinning.size()
              << " uniform enu stops defined.";
        }
        enu_range =
            std::make_pair(EnuBinning[e_stop_ctr], EnuBinning[e_stop_ctr + 1]);
      } else {
        enu_range = e_stop.get<std::pair<double, double>>("enu_range");
        if (fabs(EnuBinning.back() - enu_range.first) > 1E-6) {
          throw non_contiguous_enu_range()
              << "[ERROR]: Found enu stop {" << enu_range.first << " -- "
              << enu_range.second
              << "}, but end of last enu bin = " << EnuBinning.back()
              << ", the enu stop ranges must be contiguous.";
        }
        EnuBinning.push_back(enu_range.second);
      }
      double enu_mid = (enu_range.second + enu_range.first) / 2.0;

      std::string input_file_pattern_stop = e_stop.get<std::string>(
          "input_file_pattern",
          systtools::str_replace(input_file_pattern, "%E",
                                 StringifyNumberToOneDP(enu_mid)));
      std::string input_hist_pattern_stop = e_stop.get<std::string>(
          "input_hist_pattern",
          systtools::str_replace(input_hist_pattern, "%E",
                                 StringifyNumberToOneDP(enu_mid)));

      if (!e_stop.has_key("param_values")) {
        if (!consistent_param_values) {
          throw;
        }
        continue; // next_estop
      }

      std::vector<fhicl::ParameterSet> value_descriptors;
      size_t param_value_it = 0;
      for (fhicl::ParameterSet const &param_value :
           e_stop.get<std::vector<fhicl::ParameterSet>>("param_values")) {
        double value;
        if (!param_value.has_key("value") && consistent_param_values) {
          value = param_values[param_value_it];
        } else {
          value = param_value.get<double>("value");
        }
        param_value_it++;

        std::string input_file = param_value.get<std::string>(
            "input_file",
            systtools::str_replace(input_file_pattern_stop, "%V",
                                   StringifyNumberToOneDP(value)));
        std::string input_hist = param_value.get<std::string>(
            "input_hist",
            systtools::str_replace(input_hist_pattern_stop, "%V",
                                   StringifyNumberToOneDP(value)));

        fhicl::ParameterSet value_descriptor;
        value_descriptor.put("value", value);
        value_descriptor.put("input_file", input_file);
        value_descriptor.put("input_hist", input_hist);
        value_descriptors.push_back(std::move(value_descriptor));
      }
      estop_descriptor.put("inputs", value_descriptors);
      EnuResponses.emplace_back();
      EnuResponses.back().LoadInputHistograms(estop_descriptor);
    }
  }

  enu_bin_it_t GetEnuBin(double enu_GeV) const {
    if (EnuBinning.size() < 2) {
      return kBinOutsideRange;
    }
    if ((enu_GeV < EnuBinning.front()) || (enu_GeV >= EnuBinning.back())) {
      return kBinOutsideRange;
    }
    for (enu_bin_it_t bi_it = 0; bi_it < enu_bin_it_t(EnuBinning.size() - 1);
         ++bi_it) {
      if ((enu_GeV >= EnuBinning[bi_it]) && (enu_GeV < EnuBinning[bi_it + 1])) {
        return bi_it;
      }
    }
    return kBinOutsideRange;
  }

public:
  EnuBinnedTemplateResponseCalculator(fhicl::ParameterSet const &ps) {
    LoadInputHistograms(ps);
  };

  EnuBinnedTemplateResponseCalculator(
      EnuBinnedTemplateResponseCalculator &&other)
      : EnuBinning(std::move(other.EnuBinning)),
        EnuResponses(std::move(other.EnuResponses)) {}

  virtual std::pair<enu_bin_it_t, typename TRC::bin_it_t>
  GetBin(double enu_GeV,
         std::array<double, TRC::NDimensions> const &kinematics) const {
    enu_bin_it_t ebi_it = GetEnuBin(enu_GeV);
    if (ebi_it == kBinOutsideRange) {
      return std::pair<enu_bin_it_t, typename TRC::bin_it_t>{kBinOutsideRange,
                                                             kBinOutsideRange};
    }
    return {ebi_it, EnuResponses[ebi_it].GetBin(kinematics)};
  }
  
  // virtual destructor required
  virtual ~EnuBinnedTemplateResponseCalculator() = default;

  double
  GetVariation(double val,
               std::pair<enu_bin_it_t, typename TRC::bin_it_t> bin) const {
    return EnuResponses[bin.first].GetVariation(val, bin.second);
  }

  double
  GetVariation(double val, double enu_GeV,
               std::array<double, TRC::NDimensions> const &kinematics) const {
    return GetVariation(val, GetBin(enu_GeV, kinematics));
  }

  bool IsValidVariation(double val) {
    return EnuResponses.front().IsValidVariation(val);
  }
};

} // namespace nusyst
