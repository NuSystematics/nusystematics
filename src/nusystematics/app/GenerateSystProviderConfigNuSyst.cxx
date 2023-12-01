#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
#include "systematicstools/utility/md5.hh"

#include "nusystematics/utility/response_helper.hh"
#include "nusystematics/utility/make_instance.hh"

#include "fhiclcpp/ParameterSet.h"

#include <fstream>
#include <iomanip>
#include <iostream>

namespace cliopts {
std::string fclname = "";
std::string outputfile = "";
std::string envvar = "FHICL_FILE_PATH";
std::string fhicl_key = "syst_providers";
bool WrapWithPROLOG = true;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-c <config.fcl>  : fhicl file to read.\n"
               "\t-o <output.fcl>  : fhicl file to write, stdout by default.\n"
               "\t-k <list key>    : fhicl key to look for list of providers,\n"
               "\t                   \"syst_providers\" by default.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);
  if (!cliopts::fclname.size()) {
    std::cout << "[ERROR]: Expected to be passed a -c option." << std::endl;
    SayUsage(argv);
    exit(1);
  }

/*
  char const *ev = getenv(cliopts::envvar.c_str());
  if (!ev) {
    std::cout << "[ERROR]: Could not read environment variable:\""
              << cliopts::envvar
              << "\". Please supply a variable containing a valid path list "
                 "via the -p command line option."
              << std::endl;
    SayUsage(argv);
    exit(1);
  }

  fhicl::ParameterSet in_ps = fhicl::ParameterSet::make(cliopts::fclname, std::make_unique<cet::filepath_lookup>(ev));
*/

  // TODO
  std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();

  fhicl::ParameterSet in_ps = fhicl::ParameterSet::make(cliopts::fclname, *fm);

  std::cout << "[GenerateSystProviderConfigNuSyst] input" << std::endl;
  std::cout << in_ps.to_indented_string() << std::endl;

  std::vector<std::unique_ptr<nusyst::IGENIESystProvider_tool>> tools =
      systtools::ConfigureISystProvidersFromToolConfig<
          nusyst::IGENIESystProvider_tool>(in_ps, nusyst::make_instance,
                                           cliopts::fhicl_key);

  fhicl::ParameterSet out_ps;
  std::vector<std::string> providerNames;
  for (auto &prov : tools) {
    if (!systtools::Validate(prov->GetSystMetaData(), false)) {
      throw systtools::invalid_SystMetaData()
          << "[ERROR]: A parameter handled by provider: "
          << std::quoted(prov->GetFullyQualifiedName())
          << " failed validation.";
    }
    fhicl::ParameterSet tool_ps = prov->GetParameterHeadersDocument();
    out_ps.put(prov->GetFullyQualifiedName(), tool_ps);
    providerNames.push_back(prov->GetFullyQualifiedName());
  }
  out_ps.put("syst_providers", providerNames);

  fhicl::ParameterSet wrapped_out_ps;
  wrapped_out_ps.put("generated_systematic_provider_configuration", out_ps);

  std::ostream *os(nullptr);

  if (cliopts::outputfile.size()) {
    std::ofstream *fs = new std::ofstream(cliopts::outputfile);
    if (!fs->is_open()) {
      std::cout << "[ERROR]: Failed to open " << cliopts::outputfile
                << " for writing." << std::endl;
      exit(1);
    }
    os = fs;
  } else {
    os = &std::cout;
  }

  if (cliopts::WrapWithPROLOG) {
    (*os) << "BEGIN_PROLOG" << std::endl;
  }

  (*os) << wrapped_out_ps.to_indented_string() << std::endl;

  if (cliopts::WrapWithPROLOG) {
    (*os) << "END_PROLOG" << std::endl;
  }

  if (cliopts::outputfile.size()) {
    static_cast<std::ofstream *>(os)->close();
    delete os;
  }

  std::cout << (cliopts::outputfile.size() ? "Wrote" : "Built")
            << " systematic provider configuration with md5: "
            << std::quoted(systtools::md5(out_ps.to_compact_string())) << std::flush;
  if (cliopts::outputfile.size()) {
    std::cout << " to " << std::quoted(cliopts::outputfile) << std::flush;
  }
  std::cout << std::endl;
}
