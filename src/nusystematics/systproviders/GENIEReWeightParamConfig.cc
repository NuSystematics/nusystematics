#include "nusystematics/systproviders/GENIEReWeightParamConfig.hh"

#include "systematicstools/interface/ISystProviderTool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "systematicstools/utility/ResponselessParamUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "RwFramework/GSyst.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace systtools;
using namespace genie::rew;

namespace nusyst {

SystMetaData
ConfigureSetOfIndependentParameters(fhicl::ParameterSet const &cfg,
                                    paramId_t firstParamId,
                                    std::vector<genie::rew::GSyst_t> Dials) {

  SystMetaData MD;
  for (GSyst_t const &gdial : Dials) {
    std::string const &pname = GSyst::AsString(gdial);
    if (!FhiclToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFhiclToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    MD.push_back(std::move(param));
  }
  return MD;
}

SystMetaData ConfigureSetOfDependentShapeableParameters(
    fhicl::ParameterSet const &cfg, paramId_t firstParamId,
    fhicl::ParameterSet &tool_options, std::string const &ResponseParameterName,
    std::vector<std::pair<genie::rew::GSyst_t, genie::rew::GSyst_t>> Dials,
    bool IsShape) {

  SystMetaData MD;

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystParamHeader Response;
  if (!ignore_parameter_dependence) {
    Response.prettyName = ResponseParameterName;
    Response.systParamId = firstParamId++;
  }

  size_t NParams = 0;
  for (std::pair<GSyst_t, GSyst_t> const &gdial : Dials) {
    std::string const &pname = GSyst::AsString(gdial.first);
    if (!FhiclToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFhiclToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    param.prettyName = GSyst::AsString(IsShape ? gdial.second : gdial.first);
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = Response.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in "
               "random "
               "throw mode or set \"ignore_parameter_dependence\" "
               "in the "
               "GENIEReWeight configuration.";
      }
    }
    MD.push_back(std::move(param));
    NParams++;
  }
  if (!ignore_parameter_dependence && NParams) {
    std::vector<std::string> DialNames;
    std::transform(Dials.begin(), Dials.end(), std::back_inserter(DialNames),
                   [&](std::pair<GSyst_t, GSyst_t> const &gdial) {
                     return GSyst::AsString(IsShape ? gdial.second
                                                    : gdial.first);
                   });

    MD.push_back(std::move(Response));
    FinalizeAndValidateDependentParameters(MD, ResponseParameterName,
                                           DialNames);
  }

  return MD;
}

SystMetaData ConfigureSetOfDependentParameters(
    fhicl::ParameterSet const &cfg, paramId_t firstParamId,
    fhicl::ParameterSet &tool_options, std::string const &ResponseParameterName,
    std::vector<genie::rew::GSyst_t> Dials) {

  SystMetaData MD;

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystParamHeader Response;
  if (!ignore_parameter_dependence) {
    Response.prettyName = ResponseParameterName;
    Response.systParamId = firstParamId++;
  }

  size_t NParams = 0;
  for (GSyst_t const &gdial : Dials) {
    std::string const &pname = GSyst::AsString(gdial);
    if (!FhiclToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFhiclToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = Response.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in "
               "random "
               "throw mode or set \"ignore_parameter_dependence\" "
               "in the "
               "GENIEReWeight configuration.";
      }
    }
    MD.push_back(std::move(param));
    NParams++;
  }
  if (!ignore_parameter_dependence && NParams) {
    std::vector<std::string> DialNames;
    std::transform(Dials.begin(), Dials.end(), std::back_inserter(DialNames),
                   [](GSyst_t const &gdial) { return GSyst::AsString(gdial); });

    MD.push_back(std::move(Response));
    FinalizeAndValidateDependentParameters(MD, ResponseParameterName,
                                           DialNames);
  }

  return MD;
}

SystMetaData ConfigureQEParameterHeaders(fhicl::ParameterSet const &cfg,
                                         paramId_t firstParamId,
                                         fhicl::ParameterSet &tool_options) {
  SystMetaData QEmd;

  bool MaCCQEIsShapeOnly = cfg.get<bool>("MaCCQEIsShapeOnly", false);
  tool_options.put("MaCCQEIsShapeOnly", MaCCQEIsShapeOnly);

  // Axial FFs
  bool DipoleNormCCQEIsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormCCQE));
  bool DipoleIsShapeOnly = MaCCQEIsShapeOnly || DipoleNormCCQEIsUsed;
  bool DipoleMaCCQEIsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_MaCCQE));

  bool IsDipoleReWeight =
      DipoleIsShapeOnly || DipoleNormCCQEIsUsed || DipoleMaCCQEIsUsed;

  bool ZNormIsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZNormCCQE));
  bool ZExpA1IsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA1CCQE));
  bool ZExpA2IsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA2CCQE));
  bool ZExpA3IsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA3CCQE));
  bool ZExpA4IsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA4CCQE));

  bool IsZExpReWeight = ZNormIsUsed || ZExpA1IsUsed || ZExpA2IsUsed ||
                        ZExpA3IsUsed || ZExpA4IsUsed;

  if (IsDipoleReWeight && IsZExpReWeight) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIReWeight_tool, both dipole and "
           "Z-expansion axial form factor dials are specified.";
  }

  if (DipoleNormCCQEIsUsed && !MaCCQEIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormCCQE was "
           "requested but MaCCQE was not specified to be shape-only.";
  }

  if (IsDipoleReWeight) {
    if (DipoleNormCCQEIsUsed) {
      systtools::SystParamHeader param;
      ParseFhiclToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_NormCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
    if (DipoleMaCCQEIsUsed) {
      systtools::SystParamHeader param;
      ParseFhiclToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_MaCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      param.prettyName = GSyst::AsString(
          DipoleIsShapeOnly ? kXSecTwkDial_MaCCQEshape : kXSecTwkDial_MaCCQE);
      QEmd.push_back(std::move(param));
    }
  } else if (IsZExpReWeight) {
    // ZNorm enters in linearly to the weight calculation so doesn't need to
    // output its response via the a meta-parameter.
    if (ZNormIsUsed) {
      systtools::SystParamHeader param;
      ParseFhiclToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_ZNormCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }

    SystMetaData ZExpmd = ConfigureSetOfDependentParameters(
        cfg, firstParamId, tool_options, "ZExpAVariationResponse",
        {kXSecTwkDial_ZExpA1CCQE, kXSecTwkDial_ZExpA2CCQE,
         kXSecTwkDial_ZExpA3CCQE, kXSecTwkDial_ZExpA4CCQE});
    firstParamId += ZExpmd.size();
    ExtendSystMetaData(QEmd, std::move(ZExpmd));
  }

  // These are all independent and based upon the channel that was generated
  SystMetaData VecFFCCQEmd = ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kXSecTwkDial_VecFFCCQEshape});
  firstParamId += VecFFCCQEmd.size();
  ExtendSystMetaData(QEmd, std::move(VecFFCCQEmd));

  SystMetaData RPA_CCQEmd = ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kXSecTwkDial_RPA_CCQE});
  firstParamId += RPA_CCQEmd.size();
  ExtendSystMetaData(QEmd, std::move(RPA_CCQEmd));

  SystMetaData CoulombCCQEmd = ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kXSecTwkDial_CoulombCCQE});
  firstParamId += CoulombCCQEmd.size();
  ExtendSystMetaData(QEmd, std::move(CoulombCCQEmd));

  bool AxFFCCQEDipoleToZExp =
      cfg.get<bool>("AxFFCCQEDipoleToZExp", false) || IsZExpReWeight;

  tool_options.put<bool>("AxFFCCQEDipoleToZExp",AxFFCCQEDipoleToZExp);


  return QEmd;
}

SystMetaData ConfigureMECParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {

  return ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kXSecTwkDial_NormCCMEC,
       kXSecTwkDial_NormNCMEC,
       kXSecTwkDial_NormEMMEC,
       kXSecTwkDial_DecayAngMEC,
       kXSecTwkDial_FracPN_CCMEC,
       kXSecTwkDial_FracDelta_CCMEC,
       kXSecTwkDial_XSecShape_CCMEC}
  );

}

SystMetaData ConfigureNCELParameterHeaders(fhicl::ParameterSet const &cfg,
                                           paramId_t firstParamId,
                                           fhicl::ParameterSet &tool_options) {
  return ConfigureSetOfDependentParameters(
      cfg, firstParamId, tool_options, "NCELVariationResponse",
      {kXSecTwkDial_MaNCEL, kXSecTwkDial_EtaNCEL});
}

SystMetaData ConfigureRESParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {
  SystMetaData RESmd;

  //************* CCRES
  bool CCRESIsShapeOnly = cfg.get<bool>("CCRESIsShapeOnly", false);
  tool_options.put("CCRESIsShapeOnly", CCRESIsShapeOnly);

  bool NormCCRESIsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormCCRES));

  if (NormCCRESIsUsed && !CCRESIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormCCRES was "
           "requested but CCRES was not specified to be shape-only.";
  }

  if (NormCCRESIsUsed) {
    systtools::SystParamHeader param;
    ParseFhiclToolConfigurationParameter(
        cfg, GSyst::AsString(kXSecTwkDial_NormCCRES), param, firstParamId);
    param.systParamId = firstParamId++;
    RESmd.push_back(std::move(param));
  }

  SystMetaData CCRESmd = ConfigureSetOfDependentShapeableParameters(
      cfg, firstParamId, tool_options, "CCRESVariationResponse",
      {{kXSecTwkDial_MaCCRES, kXSecTwkDial_MaCCRESshape},
       {kXSecTwkDial_MvCCRES, kXSecTwkDial_MvCCRESshape}},
      CCRESIsShapeOnly);
  firstParamId += CCRESmd.size();
  ExtendSystMetaData(RESmd, std::move(CCRESmd));

  //************* NCRES
  bool NCRESIsShapeOnly = cfg.get<bool>("NCRESIsShapeOnly", false);
  tool_options.put("NCRESIsShapeOnly", NCRESIsShapeOnly);

  bool NormNCRESIsUsed = FhiclToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormNCRES));

  if (NormNCRESIsUsed && !NCRESIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormNCRES was "
           "requested but NCRES was not specified to be shape-only.";
  }

  if (NormNCRESIsUsed) {
    systtools::SystParamHeader param;
    ParseFhiclToolConfigurationParameter(
        cfg, GSyst::AsString(kXSecTwkDial_NormNCRES), param, firstParamId);
    param.systParamId = firstParamId++;
    RESmd.push_back(std::move(param));
  }

  SystMetaData NCRESmd = ConfigureSetOfDependentShapeableParameters(
      cfg, firstParamId, tool_options, "NCRESVariationResponse",
      {{kXSecTwkDial_MaNCRES, kXSecTwkDial_MaNCRESshape},
       {kXSecTwkDial_MvNCRES, kXSecTwkDial_MvNCRESshape}},
      NCRESIsShapeOnly);
  firstParamId += NCRESmd.size();
  ExtendSystMetaData(RESmd, std::move(NCRESmd));

  // These are all independent and based upon the channel that was generated
  SystMetaData RESOther = ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kXSecTwkDial_RvpCC1pi, kXSecTwkDial_RvpCC2pi, kXSecTwkDial_RvpNC1pi,
       kXSecTwkDial_RvpNC2pi, kXSecTwkDial_RvnCC1pi, kXSecTwkDial_RvnCC2pi,
       kXSecTwkDial_RvnNC1pi, kXSecTwkDial_RvnNC2pi, kXSecTwkDial_RvbarpCC1pi,
       kXSecTwkDial_RvbarpCC2pi, kXSecTwkDial_RvbarpNC1pi,
       kXSecTwkDial_RvbarpNC2pi, kXSecTwkDial_RvbarnCC1pi,
       kXSecTwkDial_RvbarnCC2pi, kXSecTwkDial_RvbarnNC1pi,
       kXSecTwkDial_RvbarnNC2pi,

       kRDcyTwkDial_BR1gamma, kRDcyTwkDial_BR1eta,
       kRDcyTwkDial_Theta_Delta2Npi,
       kRDcyTwkDial_Theta_Delta2NRad});
  ExtendSystMetaData(RESmd, std::move(RESOther));

  return RESmd;
}

SystMetaData ConfigureCOHParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {
  return ConfigureSetOfDependentParameters(
      cfg, firstParamId, tool_options, "COHVariationResponse",
      {kXSecTwkDial_MaCOHpi, kXSecTwkDial_R0COHpi, kXSecTwkDial_NormCCCOHpi, kXSecTwkDial_NormNCCOHpi});
}

SystMetaData ConfigureDISParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {

  bool DISBYIsShapeOnly = cfg.get<bool>("DISBYIsShapeOnly", false);
  tool_options.put("DISBYIsShapeOnly", DISBYIsShapeOnly);

  SystMetaData DISmd = ConfigureSetOfDependentShapeableParameters(
      cfg, firstParamId, tool_options, "DISBYVariationResponse",
      {{kXSecTwkDial_AhtBY, kXSecTwkDial_AhtBYshape},
       {kXSecTwkDial_BhtBY, kXSecTwkDial_BhtBYshape},
       {kXSecTwkDial_CV1uBY, kXSecTwkDial_CV1uBYshape},
       {kXSecTwkDial_CV2uBY, kXSecTwkDial_CV2uBYshape}},
      DISBYIsShapeOnly);
  firstParamId += DISmd.size();

  SystMetaData AGKYmd = ConfigureSetOfDependentParameters(
      cfg, firstParamId, tool_options, "AGKYVariationResponse",
      {kHadrAGKYTwkDial_xF1pi, kHadrAGKYTwkDial_pT1pi});
  firstParamId += AGKYmd.size();
  ExtendSystMetaData(DISmd, std::move(AGKYmd));

  SystMetaData FZmd = ConfigureSetOfIndependentParameters(
      cfg, firstParamId, {kHadrNuclTwkDial_FormZone});
  ExtendSystMetaData(DISmd, std::move(FZmd));

  return DISmd;
}

SystMetaData ConfigureFSIParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {

  SystMetaData FSImd = ConfigureSetOfDependentParameters(
      cfg, firstParamId, tool_options, "FSI_pi_VariationResponse",
      {kINukeTwkDial_MFP_pi, kINukeTwkDial_FrCEx_pi, kINukeTwkDial_FrInel_pi,
       kINukeTwkDial_FrAbs_pi, kINukeTwkDial_FrPiProd_pi});
  firstParamId += FSImd.size();

  SystMetaData FSI_N_md = ConfigureSetOfDependentParameters(
      cfg, firstParamId, tool_options, "FSI_N_VariationResponse",
      {kINukeTwkDial_MFP_N, kINukeTwkDial_FrCEx_N, kINukeTwkDial_FrInel_N,
       kINukeTwkDial_FrAbs_N, kINukeTwkDial_FrPiProd_N});

  ExtendSystMetaData(FSImd, std::move(FSI_N_md));

  return FSImd;
}

SystMetaData ConfigureOtherParameterHeaders(fhicl::ParameterSet const &cfg,
                                            paramId_t firstParamId) {
  return ConfigureSetOfIndependentParameters(
      cfg, firstParamId,
      {kSystNucl_CCQEPauliSupViaKF, kSystNucl_CCQEMomDistroFGtoSF});
}

} // namespace nusyst
