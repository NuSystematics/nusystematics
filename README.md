# `nusystematics`

Implements neutrino interaction systematics for GENIE3 events (including an 
interface to GReWeight) in the `systematicstools` framework.

## Command-line tools

A set of CLI tools under `src/nusystematics/app/` covers the typical
config-build → reweight → inspect workflow. They share a unified
`nusyst <subcommand>` dispatcher (each subcommand forwards to a dedicated
binary; both forms work):

| Subcommand | Underlying binary | Purpose |
|---|---|---|
| `nusyst config` | `GenerateAllDialsConfigNuSyst` | Kitchen-sink parameter-headers config |
| `nusyst inventory` | `DeclaredDialTestNuSyst` | Inventory + metadata dump (no event input) |
| `nusyst tweaks` | `DumpConfiguredTweaksNuSyst` | Flat TTree of per-event weights (`-j N` for multi-process) |
| `nusyst response` | `DumpDialResponseNuSyst` | Weight-vs-paramValue curves per dial |
| `nusyst plots` | `PlotSystVariationsNuSyst` | Differential-xsec plots with ratio panels, split by channel |

`nusyst help` lists the subcommands; `nusyst <subcmd> -?` shows tool help.

See [`doc/CLI_TOOLS.md`](doc/CLI_TOOLS.md) for environment setup, the
end-to-end recipe, per-tool reference, and the `applies_to_channels` /
channel-bucket architecture used to skip dial × channel combinations the
dial doesn't apply to.

## To Build

`nusystematics` will by default build `systematicstools` for you, but requires 
ROOT v6+ and GENIE v3+ to be set up in the environment at the time of 
configuration. Requires cmake 3.21+ (see [here](#fetch-cmake-321) for updating CMake hints).

Having checked out the repository in `/path/to/repo`:

```bash
cd /path/to/repo
mkdir build; cd build
cmake ..
make install -j 8
```

## Adding a new systprovider

First off, read documentation describing systematicstools design choices 
[here](https://github.com/jedori0228/systematicstools#introduction).

Configuration is 'two step', first from a [_tool configuration_](https://github.com/jedori0228/systematicstools/blob/develop/doc/ToolConfiguration.md#tool-configuration) 
FHiCL file, which produces the [_parameter headers_](https://github.com/jedori0228/systematicstools/blob/develop/doc/ParameterHeaders.md) 
FHiCL file that is used for the second step. The parameter headers file is then 
used to configure an instance of your systprovider into a state where it can 
produce event responses. For more high-level information on writing a provider, 
see [here](https://github.com/jedori0228/systematicstools/blob/develop/doc/WritingAProvider.md).

## Learning the package

See [DUNE Data Analysis School 2026 matrial](https://github.com/jedori0228/DDAS2026_NuSyst).