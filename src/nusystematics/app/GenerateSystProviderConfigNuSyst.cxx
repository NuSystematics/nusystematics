#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/make_instance.hh"

#include "yaml-cpp/yaml.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace cliopts {
std::string yamlname = "";
std::string outputfile = "";
std::string envvar = "YAML_FILE_PATH";
std::string yaml_key = "syst_providers";
bool WrapWithPROLOG = false;
} // namespace cliopts

namespace {

void SetAllSequencesToFlow(YAML::Node node) {
  if (!node) {
    return;
  }

  if (node.IsSequence()) {
    node.SetStyle(YAML::EmitterStyle::Flow);
    for (auto const &child : node) {
      SetAllSequencesToFlow(child);
    }
    return;
  }

  if (node.IsMap()) {
    for (auto const &kv : node) {
      SetAllSequencesToFlow(kv.second);
    }
  }
}

// Recursively load a YAML file, merging any files listed under an
// "includes" key (paths resolved relative to the including file).
// Keys from the including file overwrite keys from included files.
// The "includes" key itself is consumed and not forwarded.
YAML::Node LoadYAMLWithIncludes(std::string const &filepath) {
  YAML::Node node = YAML::LoadFile(filepath);

  if (!node["includes"]) {
    return node;
  }

  // base_dir: the directory containing the current file, used to resolve
  // relative paths in "includes". E.g. if filepath is "/a/b/top.yaml",
  // base_dir is "/a/b", so "sub/foo.yaml" resolves to "/a/b/sub/foo.yaml".
  std::filesystem::path base_dir =
      std::filesystem::path(filepath).parent_path();

  // Start with merged content from all included files (in order).
  YAML::Node merged;
  for (auto const &inc : node["includes"]) {
    // inc_path: path to the included file as written in the YAML.
    // Environment variables in ${VAR} syntax are expanded first, then if the
    // result is still a relative path it is resolved against base_dir so the
    // lookup is always relative to the including file, not the working directory.
    std::filesystem::path inc_path(systtools::expand_env_vars(inc.as<std::string>()));
    if (inc_path.is_relative()) {
      inc_path = base_dir / inc_path;
    }
    YAML::Node inc_node = LoadYAMLWithIncludes(inc_path.string());
    for (auto const &kv : inc_node) {
      std::string key = kv.first.as<std::string>();
      if (key == "includes") continue;
      merged[key] = kv.second;
    }
  }

  // Overlay the top-level file's own keys (except "includes").
  // This lets the top-level file override anything from included files,
  // including "syst_providers".
  for (auto const &kv : node) {
    std::string key = kv.first.as<std::string>();
    if (key == "includes") continue;
    merged[key] = kv.second;
  }

  return merged;
}

} // namespace

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-c <config.yaml> : YAML file to read.\n"
               "\t-o <output.yaml> : YAML file to write, stdout by default.\n"
               "\t-k <list key>    : YAML key to look for list of providers,\n"
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
      cliopts::yamlname = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::yaml_key = argv[++opt];
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
  if (!cliopts::yamlname.size()) {
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

    YAML::Node in_ps = YAML::Node::make(cliopts::yamlname,
    std::make_unique<cet::filepath_lookup>(ev));
  */

  YAML::Node in_yaml = LoadYAMLWithIncludes(cliopts::yamlname);

  std::cout << "[GenerateSystProviderConfigNuSyst] input" << std::endl;
  std::cout << in_yaml << std::endl;

  std::vector<std::unique_ptr<nusyst::IGENIESystProvider_tool>> tools =
      systtools::ConfigureISystProvidersFromToolConfig<
          nusyst::IGENIESystProvider_tool>(in_yaml, nusyst::make_instance,
                                           cliopts::yaml_key);

  YAML::Node out_yaml;
  std::vector<std::string> providerNames;
  for (auto &prov : tools) {
    if (!systtools::Validate(prov->GetSystMetaData(), false)) {
      throw systtools::invalid_SystMetaData()
          << "[ERROR]: A parameter handled by provider: "
          << std::quoted(prov->GetFullyQualifiedName())
          << " failed validation.";
    }
    YAML::Node tool_yaml = prov->GetParameterHeadersDocument();
    out_yaml[prov->GetFullyQualifiedName()] = tool_yaml;
    providerNames.push_back(prov->GetFullyQualifiedName());
  }
  out_yaml["syst_providers"] = providerNames;
  out_yaml["syst_providers"].SetStyle(YAML::EmitterStyle::Flow);

  YAML::Node wrapped_out_yaml;
  wrapped_out_yaml["generated_systematic_provider_configuration"] = out_yaml;
  SetAllSequencesToFlow(wrapped_out_yaml);

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
    (*os) << "# YAML Configuration" << std::endl;
  }

  (*os) << wrapped_out_yaml << std::endl;

  if (cliopts::outputfile.size()) {
    static_cast<std::ofstream *>(os)->close();
    delete os;
  }

  std::string yaml_str = YAML::Dump(out_yaml);
  std::cout << (cliopts::outputfile.size() ? "Wrote" : "Built")
            << " systematic provider configuration with md5: "
            << std::quoted(systtools::md5(yaml_str))
            << std::flush;
  if (cliopts::outputfile.size()) {
    std::cout << " to " << std::quoted(cliopts::outputfile) << std::flush;
  }
  std::cout << std::endl;
}
