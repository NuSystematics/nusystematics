#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/make_instance.hh"
#include "nusystematics/utility/silence_genie.hh"

#include "yaml-cpp/yaml.h"

#include <cstdlib>
#include <filesystem>
#include "Framework/Messenger/Messenger.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace cliopts {
std::string yamlname = "";
std::string outputfile = "";
std::string envvar = "YAML_FILE_PATH";
std::string yaml_key = "syst_providers";
bool WrapWithPROLOG = false;
bool DoDebug = false;
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

std::vector<std::filesystem::path> GetYamlSearchPaths() {
  std::vector<std::filesystem::path> search_paths;

  char const *config_dirs = std::getenv("NUSYST_CONFIG_DIR");
  if (!config_dirs || !*config_dirs) {
    return search_paths;
  }

  std::stringstream dirs_stream(config_dirs);
  std::string dir;
  while (std::getline(dirs_stream, dir, ':')) {
    if (!dir.empty()) {
      search_paths.emplace_back(dir);
    }
  }

  return search_paths;
}

std::filesystem::path ResolveYamlPath(std::string const &filepath,
                                      std::filesystem::path const &base_dir = {}) {
  std::filesystem::path path(systtools::expand_env_vars(filepath));

  auto const exists = [](std::filesystem::path const &candidate) {
    return !candidate.empty() && std::filesystem::exists(candidate);
  };

  if (path.is_absolute()) {
    if (exists(path)) {
      return path;
    }
  } else {
    if (!base_dir.empty()) {
      std::filesystem::path candidate = base_dir / path;
      if (exists(candidate)) {
        return candidate;
      }
    }

    if (exists(path)) {
      return path;
    }

    for (auto const &search_dir : GetYamlSearchPaths()) {
      std::filesystem::path candidate = search_dir / path;
      if (exists(candidate)) {
        return candidate;
      }
    }
  }

  throw std::runtime_error("Unable to locate YAML file: " + path.string());
}

// Recursively load a YAML file, merging any files listed under an
// "includes" key (paths resolved relative to the including file).
// Keys from the including file overwrite keys from included files.
// The "includes" key itself is consumed and not forwarded.
YAML::Node LoadYAMLWithIncludes(std::string const &filepath,
                                std::filesystem::path const &base_dir = {}) {
  std::filesystem::path resolved_filepath = ResolveYamlPath(filepath, base_dir);
  YAML::Node node = YAML::LoadFile(resolved_filepath.string());

  if (!node["includes"]) {
    return node;
  }

  // base_dir: the directory containing the current file, used to resolve
  // relative paths in "includes". E.g. if filepath is "/a/b/top.yaml",
  // base_dir is "/a/b", so "sub/foo.yaml" resolves to "/a/b/sub/foo.yaml".
  std::filesystem::path next_base_dir = resolved_filepath.parent_path();

  // Start with merged content from all included files (in order).
  YAML::Node merged;
  for (auto const &inc : node["includes"]) {
    // inc_path: path to the included file as written in the YAML.
    // Environment variables in ${VAR} syntax are expanded first, then if the
    // result is still a relative path we try the including file's directory
    // first, then the directories listed in NUSYST_CONFIG_DIR.
    YAML::Node inc_node = LoadYAMLWithIncludes(inc.as<std::string>(), next_base_dir);
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
               "\t--debug         : Run debug mode.\n"
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
    } else if (std::string(argv[opt]) == "--debug") {
      cliopts::DoDebug = true;
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

  YAML::Node in_yaml = LoadYAMLWithIncludes(cliopts::yamlname);

  if(cliopts::DoDebug){
    std::cout << "[GenerateSystProviderConfigNuSyst] input" << std::endl;
    std::cout << in_yaml << std::endl;
  }

  nusyst::quiet::SetGlobalQuiet();

  std::vector<std::unique_ptr<nusyst::IGENIESystProvider_tool>> tools;
  {
    // Suppress GENIE/provider chatter during config load and provider build.
    nusyst::quiet::StdoutSink _quiet;
    genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
        "Messenger_whisper.xml");

    tools = systtools::ConfigureISystProvidersFromToolConfig<
            nusyst::IGENIESystProvider_tool>(in_yaml, nusyst::make_instance,
                                           cliopts::yaml_key);
  }

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
