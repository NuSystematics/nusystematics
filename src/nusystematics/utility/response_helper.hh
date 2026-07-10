#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/GENIEReWeight_tool.hh"

#include "nusystematics/utility/make_instance.hh"

#include "systematicstools/interface/SystParamHeader.hh"

#include "systematicstools/interpreters/ParamHeaderHelper.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"

#include "Framework/EventGen/EventRecord.h"

#include "nusystematics/utility/Flatness.h" //header added ~MM

#include "nusystematics/utility/RoundRobin.h" //header added ~MM

#include "TFile.h"
#include "TTree.h"

#include <chrono>
#include <string>

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(response_helper_found_no_parameters);

class response_helper : public systtools::ParamHeaderHelper {

  size_t NEvsProcessed;
  std::map<simb_mode_copy, std::map<size_t, std::tuple<double, double, size_t>>>
      ProfileStats;

private:
  constexpr static size_t Order = 5;
  constexpr static size_t NCoeffs = Order + 1;

  size_t ProfilerRate;

  std::string config_file;
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> syst_providers;

public:
  response_helper() : NEvsProcessed(0), ProfilerRate(0) {}
  response_helper(std::string const &fhicl_config_filename) : NEvsProcessed(0) {
    LoadConfiguration(fhicl_config_filename);
  }

  std::vector<std::unique_ptr<IGENIESystProvider_tool>>& GetSystProvider(){
    return syst_providers;
  };

  void LoadProvidersAndHeaders(fhicl::ParameterSet const &ps) {
    syst_providers = systtools::ConfigureISystProvidersFromParameterHeaders<
        IGENIESystProvider_tool>(ps, make_instance);
    
    if (!syst_providers.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected to load some systematic providers from input: "
          << std::quoted(config_file);
    }
    
    systtools::param_header_map_t configuredParameterHeaders =
        systtools::BuildParameterHeaders(syst_providers);
    if (!configuredParameterHeaders.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected systematric providers loaded from input: "
          << std::quoted(config_file) << " to provide some parameter headers.";
    }
    
    SetHeaders(configuredParameterHeaders);
  }

  void LoadConfiguration(std::string const &fhicl_config_filename) {
    config_file = fhicl_config_filename;

    // TODO
    std::unique_ptr<cet::filepath_maker> fm
        = std::make_unique<cet::filepath_lookup_nonabsolute>("FHICL_FILE_PATH");
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(config_file, *fm);

    LoadProvidersAndHeaders(ps.get<fhicl::ParameterSet>(
        "generated_systematic_provider_configuration"));

    ProfilerRate = ps.get<size_t>("ProfileRate", 0);
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord const &GenieGHep,
                    systtools::paramId_t i) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep, i);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  systtools::event_unit_response_w_cv_t
  GetEventVariationAndCVResponse(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_w_cv_t response;

    simb_mode_copy mode = GetSimbMode(GenieGHep);

    for (size_t sp_it = 0; sp_it < syst_providers.size(); ++sp_it) {
      std::unique_ptr<IGENIESystProvider_tool> const &sp =
          syst_providers[sp_it];

      std::chrono::high_resolution_clock::time_point start;
      if (ProfilerRate) {
        start = std::chrono::high_resolution_clock::now();
      }

      systtools::event_unit_response_w_cv_t prov_response =
          sp->GetEventVariationAndCVResponse(GenieGHep);

      if (ProfilerRate && prov_response.size()) {
        auto end = std::chrono::high_resolution_clock::now();
        auto diff_ms =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        if (!ProfileStats[mode].count(sp_it)) {
          ProfileStats[mode][sp_it] = {0, 0, 0};
        }

        std::get<2>(ProfileStats[mode][sp_it])++;

        double delta = diff_ms - std::get<0>(ProfileStats[mode][sp_it]);

        std::get<0>(ProfileStats[mode][sp_it]) +=
            delta / std::get<2>(ProfileStats[mode][sp_it]);

        std::get<1>(ProfileStats[mode][sp_it]) +=
            delta * (diff_ms - std::get<0>(ProfileStats[mode][sp_it]));
      }
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }

    if (ProfilerRate && NEvsProcessed && !(NEvsProcessed % ProfilerRate)) {
      std::cout << std::endl
                << "[PROFILE]: Event number = " << NEvsProcessed << std::endl;
      for (size_t sp_it = 0; sp_it < syst_providers.size(); ++sp_it) {
        std::cout << "\tSystProvider: "
                  << syst_providers[sp_it]->GetFullyQualifiedName()
                  << std::endl;
        for (auto const &m : ProfileStats) {
          if (!m.second.count(sp_it) || (std::get<2>(m.second.at(sp_it)) < 2)) {
            continue;
          }
          std::cout << "\t\tMode: " << tostr(m.first)
                    << ", NEvs = " << std::get<2>(m.second.at(sp_it))
                    << ", mean: " << std::get<0>(m.second.at(sp_it))
                    << " us, stddev: "
                    << sqrt(std::get<1>(m.second.at(sp_it)) /
                            (std::get<2>(m.second.at(sp_it)) - 1))
                    << " us." << std::endl;
        }
      }
    }
    NEvsProcessed++;

    return response;
  }

  double GetEventWeightResponse(genie::EventRecord const &GenieGHep,
                                systtools::param_value_list_t const &vals) {
    double weight = 1;
    for (auto &sp : syst_providers) {
      weight *= sp->GetEventWeightResponse(GenieGHep, vals);
    }
    return weight;
  }

  // Improved spline-based response
  std::vector<std::vector<double>> GetImprovedParameterResponse(
    systtools::paramId_t pid, 
    genie::EventRecord const &GenieGHep){

    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] syst_providers.size() = " << syst_providers.size() << std::endl;
    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] pid = " << pid << std::endl;

    // Get SystProvider that contains param with this given pid
    size_t sp_idx = systtools::kParamUnhandled<size_t>;
    for(size_t i=0; i<syst_providers.size(); i++){
      if(syst_providers[i]->ParamIsHandled(pid)){
        sp_idx = i;
        break;
      }
    }
    if(sp_idx==systtools::kParamUnhandled<size_t>){
      throw response_helper_found_no_parameters()
          << "[ERROR]: ParamID = " << pid << " is not found from configured parameter headers";
    }

    std::unique_ptr<IGENIESystProvider_tool> const &sp = syst_providers[sp_idx];
    IGENIESystProvider_tool* rawPtr = sp.get();
    GENIEReWeight* genie_sp = dynamic_cast<GENIEReWeight*>(rawPtr);

    // Now this only works for GENIEReWeight tool
    // if(sp->GetToolType()!="GENIEReWeight"){
    //   return GetParameterResponse(pid, v, eur);
    // }

    // Getting the header
    std::vector<double> og_nodes;
    systtools::SystParamHeader const & sph = GetHeader(pid);
    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] Prining variations set by fcl" << std::endl;
    for(const auto& v: sph.paramVariations){
      //std::cout << v << std::endl;
      og_nodes.push_back(v);
    }

    // now we know i-th systprovider is GENIEReWeight
    // get the GENIEResponseParameter with this given pid
    std::vector<std::string> rw_calc_names;
    GENIEResponseParameter& genieRP = genie_sp->GetGENIEResponseParameter(pid);
    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] Number of GReWeights: " << genieRP.Herg.size() << std::endl;
    for(auto& grw: genieRP.Herg){
    // grw is our std::unique_ptr<genie::rew::GReWeight>
    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] - Printing Infos of this GReWeight" << std::endl;
    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse]   - WghtCalcNames" << std::endl;
    for(auto& name: grw->WghtCalcNames()){
      //std::cout << name << std::endl;
    }
      //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse]   - GSystSet" << std::endl;
    genie::rew::GSystSet &gss = grw->Systematics();
    for(auto& gs: gss.AllIncluded()){
    //     std::string gs_string = genie::rew::GSyst::AsString(gs);
    //     //printf("%s %d\n", gs_string.c_str(), gs);
     }

    }

    //std::cout << "[JSKIMDEBUG][GetImprovedParameterResponse] Printing GSysts" << std::endl;
    for(auto& grw: genieRP.Herg){
      // grw is our std::unique_ptr<genie::rew::GReWeight
      for(auto& name: grw->WghtCalcNames()){
        //std::cout << name << std::endl;
      }
    }

    //BEGIN TEST AND FIX
    bool passed=false;
    int iteration=0;
    std::vector<std::vector<double>> response_vec;
    //This is where initial parameter node points and weights are stored
    //Ensures that the whatever is stored is not altered, so the vectors below will be alterd
    std::vector<double> og_weights; //= eur[0].responses;
    //Scrap above and use the below code to set up the weights since the method below returns correct values
    //std::cout << "[MazenMalak] size of file of systematics " << genieRP.Herg.size() << std::endl;
    for(int i=0; i<og_nodes.size(); i++){
      genieRP.Herg[0]->Systematics().Set(genieRP.Herg[0]->Systematics().AllIncluded()[0],og_nodes[i]);
      genieRP.Herg[0]->Reconfigure();
      og_weights.push_back(genieRP.Herg[0]->CalcWeight(GenieGHep));
    }
    // genieRP.Herg[pid]->Systematics().Set(genieRP.Herg[pid]->Systematics().AllIncluded()[pid],-3.0);
    // genieRP.Herg[pid]->Reconfigure();
    // std::cout << "ASK GENIE (1)" <<  genieRP.Herg[pid]->CalcWeight(GenieGHep) << std::endl;
    // for(int i=0; i<eur[0].responses.size(); i++){
    //   // og_weights.push_back(eur[0].responses[i]);
    //   std::cout << "[Mazen_Malak resp by GENIE] node, wieght--> " << og_nodes[i] << ", " << og_weights[i] << std::endl;
    //   std::cout << "ASK GENIE (2)" <<  genieRP.Herg[pid]->CalcWeight(GenieGHep) << std::endl;
    // }
    // std::cout << "ASK GENIE (3)" <<  genieRP.Herg[pid]->CalcWeight(GenieGHep) << std::endl;
    // std::cout << "Improved SPlining Alg" << std::endl;
    while(passed==false && iteration<10){
      if(iteration>0){ //we dont want to fix for fully flat splines (causes alot of trouble so put first iter through round robin)
        std::vector<double> Piecewise_check = Flatness(GenieGHep, pid, genieRP.Herg[0],og_nodes,og_weights);
        //if flatness check returns a size of 1, it was never flat
        if(Piecewise_check.size()!=1){
          //proabbaly a better way to do this, but this is the erasing and addtion
          //method of new nodes/weights to the vector
          int genie_stop_index=0;
          for(int i=0; i<og_nodes.size(); i++){
              if(og_nodes[i]<Piecewise_check[0]){
                  genie_stop_index++;
              }
          }
          //std::cout << "flat event" << std::endl;
          //std::cout << Piecewise_check[0] << std::endl; 
          if(Piecewise_check[0]>0){ //flatness on right side
            og_nodes.erase(og_nodes.begin()+genie_stop_index, og_nodes.end());
            og_weights.erase(og_weights.begin()+genie_stop_index, og_weights.end());
            og_nodes.push_back(Piecewise_check[0]);
            og_weights.push_back(Piecewise_check[1]);
          }
          else{ //flatness on left side
            og_nodes.erase(og_nodes.begin(), og_nodes.begin()+genie_stop_index);
            og_weights.erase(og_weights.begin(), og_weights.begin()+genie_stop_index);
            og_nodes.insert(og_nodes.begin(), Piecewise_check[0]);
            og_weights.insert(og_weights.begin(), Piecewise_check[1]);
          }
          //std::cout << Piecewise_check[1] << std::endl;
        }
      }
      //begin the round robin test
      std::array<double, 3> RoundRobin_Test = RoundRobin(og_nodes, og_weights);
      if(RoundRobin_Test[0]==0.0){
        //Find index of acquired max point deviation and append the node points to the left
        //and right of the max point along with the asociated weight]
        auto it = find(og_nodes.begin(), og_nodes.end(), RoundRobin_Test[1]);
        int index = it-og_nodes.begin();
        double left_point=RoundRobin_Test[1]-0.5*abs(RoundRobin_Test[1]-og_nodes[index-1]);
        double right_point=RoundRobin_Test[1]+0.5*abs(RoundRobin_Test[1]-og_nodes[index+1]);
        og_nodes.push_back(left_point);
        og_nodes.push_back(right_point);
        sort(og_nodes.begin(), og_nodes.end());
        if(og_nodes[0]==og_nodes[1]){
            og_nodes.erase(og_nodes.begin());
        }
        auto it_r = find(og_nodes.begin(), og_nodes.end(), right_point);
        auto it_l = find(og_nodes.begin(), og_nodes.end(), left_point);
        int index_r = it_r-og_nodes.begin();
        int index_l = it_l-og_nodes.begin();
        genieRP.Herg[0]->Systematics().Set(genieRP.Herg[0]->Systematics().AllIncluded()[0],right_point);//.Set(kINukeTwkDial_FrAbs_pi, right_point); //Change per syst
        genieRP.Herg[0]->Reconfigure();
        og_weights.insert(og_weights.begin()+index_r-1, genieRP.Herg[0]->CalcWeight(GenieGHep));
        genieRP.Herg[0]->Systematics().Set(genieRP.Herg[0]->Systematics().AllIncluded()[0], left_point); //Change per syst
        genieRP.Herg[0]->Reconfigure();
        og_weights.insert(og_weights.begin()+index_l, genieRP.Herg[0]->CalcWeight(GenieGHep));
        iteration++;
        //std::cout << "max dev= " << RoundRobin_Test[2] << " : new point added at param =" << right_point << " and " << left_point << std::endl;
      }
      else{
        passed=true;
        response_vec.push_back(og_nodes);
        response_vec.push_back(og_weights);
        //std::cout << response_vec[0].size() << " by " << response_vec[1].size() << std::endl;
      }
      if(iteration>9 && passed==false){
        response_vec.push_back(og_nodes);
        response_vec.push_back(og_weights);
      }
    }
    return response_vec;
  }

}; // class def response_helper 

} // namespace nusyst
