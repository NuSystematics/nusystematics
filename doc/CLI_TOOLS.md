# nusystematics CLI tools

This document describes the command-line tools shipped under
`src/nusystematics/app/`. They cover four stages of the typical workflow:

| Stage | Subcommand | Underlying binary | Purpose |
|---|---|---|---|
| 1. Build a config | `nusyst config` | `GenerateAllDialsConfigNuSyst` | Emit a parameter-headers fhicl with every available dial |
| 2. Inspect | `nusyst inventory` | `DeclaredDialTestNuSyst` | Print declared dials + metadata; no event input |
| 3. Reweight events | `nusyst tweaks` | `DumpConfiguredTweaksNuSyst` | Run every provider on a GHEP ntuple and write a flat TTree of responses |
| 4. Visualise | `nusyst plots` | `PlotSystVariationsNuSyst` | Per-channel differential-xsec plots with ratio panels |
| 4'. Visualise | `nusyst response` | `DumpDialResponseNuSyst` | Per-channel weight-vs-paramValue curves for every dial |

`nusyst` is a thin shell dispatcher installed alongside the binaries. Each
subcommand forwards its arguments verbatim to the underlying tool — they
remain callable by their full name if you prefer. `nusyst help` lists the
subcommands; `nusyst <subcmd> -?` shows the tool-specific help.

All binaries link `nusyst::systproviders` and `nusyst::utility` and depend
on a working GENIE 3+ / GENIE-Reweight / ROOT 6+ environment at runtime.
Source `<install-prefix>/bin/setup.nusystematics.sh` after `make install` to
put the binaries (and `nusyst`) on `PATH` and `$NUSYST/fcl` on `FHICL_FILE_PATH`.

---

## Kitchen-sink end-to-end recipe

To just see *what dials are available* on this build, with no setup:

```bash
nusyst inventory
```

The tool auto-generates a kitchen-sink config (`--mode all`) into `/tmp`,
runs the inventory, and deletes the temp file on exit. Pass `--keep-tmp` if
you want to inspect the generated fhicl. Pass `-c <file.fcl>` to inventory
an existing config instead.

For the full reweight workflow against a real GHEP ntuple:

```bash
GHEP=<path/to/ghep.root>

# 1. Build a kitchen-sink dial config:
#    --mode all = every GENIE Reweight dial (split by channel bucket) PLUS
#                 every standalone provider tool config under $NUSYST/fcl
#                 that loads successfully — incl. ZExpPCAWeighter's b1..b4
#                 (covariance-rotated CCQE axial-FF dials, see "Z-expansion
#                 dials" below) and ResIso.
nusyst config    --mode all -o /tmp/all.fcl

# 2. Run every dial on every event into a flat tree (multi-process: 4 workers).
nusyst tweaks    -c /tmp/all.fcl -i $GHEP \
                 -o /tmp/all_dump.root -N 5000 -j 4

# 3a. Per-event weight-vs-paramVal curves split by channel.
nusyst response  -c /tmp/all.fcl -i $GHEP \
                 -N 5000 -n 5 -o /tmp/all_dialresp

# 3b. Differential-xsec plots with ratio panels, split by channel.
nusyst plots     -c /tmp/all.fcl -i /tmp/all_dump.root \
                 -o /tmp/all_plot
```

> **Heads up — both CCQE Z-expansion dial sets are present in `--mode all`.**
> The raw `ZExpA1..4CCQE` + `ZNormCCQE` come from `GENIEReWeight_CCQE` and are
> *uncorrelated* — each shifts one coefficient by its own per-coefficient
> fractional error. The PCA dials `b1..b4` come from `ZExpPCAWeighter` and
> rotate into the eigenbasis of the published PRD 93, 113015 covariance
> matrix, so they're truly independent unit-variance coordinates. Both are
> exposed for validation/inspection purposes — **do not enable both groups
> simultaneously in a fit or systematic envelope**: they parameterise the
> same physics in two different bases and would double-count.

`PlotSystVariationsNuSyst` also accepts a raw GHEP + tool fhicl directly,
skipping step 3 — input format is auto-detected.

---

## Tool reference

### `GenerateAllDialsConfigNuSyst`

Builds a parameter-headers fhicl that declares every dial reachable from the
current install. Three modes:

| `--mode` | Source of dials |
|---|---|
| `genierw` | Every `genie::rew::GSyst_t` enumerated at runtime via `GSyst::AsString` |
| `providers` | Every standalone provider tool config under `--fcl-dir` (default `$NUSYST/fcl`) |
| `all` *(default)* | Both |

Key flags:

- `--single-instance` — emit one big `GENIEReWeight_All` provider with every
  dial. The default splits dials across per-channel buckets
  (`CCQE`, `CCRES`, `NCEL`, `NCRES`, `RES`, `COH`, `DIS`, `MEC`, `SPP`, `FSI`).
  See [Channel buckets and `applies_to_channels`](#channel-buckets-and-applies_to_channels)
  for what the split unlocks downstream.
- `--include-skeleton` — keep `SkeleWeighter.ToolConfig.fcl` and other
  development templates that are skipped by default. SkeleWeighter currently
  segfaults due to an upstream GReWeight rule-of-three violation; only use
  this if you are actively developing SkeleWeighter.
- `--skip <a,b,...>` — comma-separated tool-config filenames to skip in
  `providers` mode (on top of the default skip list).
- `--variation-descriptor "[-3,...,3]"` — variation grid used for every
  GENIE RW dial. Default `[-3, -2, -1, 0, 1, 2, 3]`.
- `-o <out.fcl>` — output file (default stdout).

#### Dials excluded by default

The default `genierw` mode prunes several GENIE Reweight dials before writing.
Run the tool once and look at the `Skipped GENIE RW dials:` section for the
current list and per-dial reason:

- **The wrong CCQE axial-FF dial family for the loaded tune** — GENIE's CCQE
  axial form factor has two parameterisations: a 1-parameter dipole form
  knobbed by `MaCCQE` / `E0CCQE`, and a multi-parameter Z-expansion form
  knobbed by `ZNormCCQE` + `ZExpA1..4CCQE`. A `GReWeightNuXSecCCQE` engine
  auto-selects its mode at construction by reading
  `FormFactorsAlg/AxialFormFactorModel` from the loaded tune's CCQE algorithm
  config — dipole tunes get `kModeNormAndMaShape` / `kModeMa`, Z-expansion
  tunes (e.g. AR23) get `kModeZExp`. The engine then accepts only the
  matching dial family — the wrong-mode dial gets a
  `WARN ReW: Systematic X is not handled for algorithm Y` and is silently
  ignored. `GenerateAllDialsConfigNuSyst` probes the engine the same way at
  startup, prints the detected mode in a `[INFO] CCQE axial form factor
  detected: ...` line, and drops the inactive family from the output config
  with reason `"dipole CCQE FF — loaded tune uses Z-expansion axial FF"`
  (or its mirror). If `GENIE_XSEC_TUNE` isn't set the probe returns
  `unknown` and both families are emitted.
  See also the **`ZExpPCAWeighter` (PCA-rotated b₁..b₄)** discussion below —
  on a Z-expansion tune that provider is the covariance-correct way to vary
  the same FF, and a `--mode all` config will declare both the raw
  `ZExpA*CCQE` and the PCA `b*` dials so the user can compare; don't enable
  both groups simultaneously in a fit.
- **Shape-only twins** (`MaCCQEshape`, `Ma*RESshape`, `Mv*RESshape`,
  `AhtBYshape`, `BhtBYshape`, `CV1uBYshape`, `CV2uBYshape`, `E0CCQEshape`,
  `VecFFCCQEshape`) — redundant with the shape+norm versions and unsafe on
  non-matching topologies.
- **`Norm{CCQE,CCRES,NCRES}`** — require the corresponding `IsShapeOnly`
  tool option, which conflicts with keeping the shape+norm dials.
- **EM-only dials** (`NormEMMEC`, `EmpMEC_FracPN_EM`, `EmpMEC_FracEMQE`) —
  only fire on EM probe scattering; trivial on every weak-current neutrino event.

#### `[FAIL]` vs `[SKIP]` in `providers` mode

When scanning `--fcl-dir` the tool tries to instantiate each tool config. Two
outcomes are surfaced:

- `[SKIP]` — instantiation threw an error that looks like a missing-data
  portability issue (e.g. an absolute path baked into the tool config that
  doesn't exist on this host). The dial isn't a code bug — it just needs the
  data file in place. Reported but doesn't fail the run.
- `[FAIL]` — instantiation threw for any other reason. The dial is genuinely
  broken in this build and is excluded from the output config.

---

### `DumpConfiguredTweaksNuSyst`

Reads a parameter-headers fhicl + a GHEP TChain and writes a flat TTree
(`events`) of per-event responses for every configured dial, plus a small
`tweak_metadata` TTree describing the dial names and tweak values.

```bash
DumpConfiguredTweaksNuSyst -c <config.fcl> -i <ghep.root> -o <out.root> \
                            [-N <NMax>] [-s <NSkip>] [-j N]
```

Key flags:

- `-j N` / `--threads N` — fork `N` worker processes, partition the input
  events into `N` disjoint ranges, and `hadd` the parts. Defaults to 1 (serial).
  GENIE's reweight stack holds non-thread-safe globals (`Messenger`, `RandomGen`),
  so the implementation uses `fork()`, not `std::thread`. On a 1000-event
  AR23 ν<sub>μ</sub>+Ar40 sample, `-j 4` is ~2.4× faster than `-j 1`.
- `-k <key>` — top-level fhicl key (default
  `generated_systematic_provider_configuration`).
- `-b <branch>` — `NtpMCEventRecord` branch name (default `gmcrec`).

#### Channel-aware skip

If the fhicl was produced by `GenerateAllDialsConfigNuSyst` *without*
`--single-instance`, each provider block carries a
`tool_options.applies_to_channels` list of fnmatch glob patterns (see
[Channel buckets and `applies_to_channels`](#channel-buckets-and-applies_to_channels)).
For every event, providers whose patterns don't match the event's channel key
are skipped entirely — they get a synthesised trivial response (CV=1, every
variation=1.0) so the output schema stays stable. Reported on startup:

```
[INFO]: 10 of 11 providers declare applies_to_channels;
        1 will be evaluated on every event.
```

The single un-filtered provider in that line is `FSI` — FSI reweighters
apply to any topology with FS hadrons traversing the nucleus.

#### GENIE-version detection

At startup the tool reads the `genie::NtpMCTreeHeader` record from the input
file and prints the `cvstag` (e.g. `R-3_04_00`) and the tune string. If the
sample uses a "999.999.999" placeholder, the line reads
`GENIE sample version: could not detect`. Useful when chasing version-skew
errors like the Jacobian crashes that hit GENIE 3.04 QE events under 3.06+
reweighters.

---

### `PlotSystVariationsNuSyst`

Produces per-channel differential-xsec plots with ratio panels, plus
inclusive Total / Total_CC / Total_NC pages and a one-page summary of event
counts. Accepts either a flat-tree dump from `DumpConfiguredTweaksNuSyst`
*or* a raw GHEP file directly (input format is auto-detected).

```bash
PlotSystVariationsNuSyst -i <input.root> [-c <config.fcl>] [-o <out>] \
                          [-v Enu,Q2,...] [-p MaCCQE,FormZone] [-N <NMax>] \
                          [--nbins N] [--no-pdf | --no-root]
```

Predefined kinematic variables (used by default if `-v` isn't given):
`Enu, q0, Q2, q3, W, plep, coslep, EAvail, Tp, Emiss, pmiss, Ereco_cal,
nproton, npip, npim, npi0, nneutron`. Variables whose branches aren't in the
input are silently skipped with a `[WARN]:` line.

`-c` is **required for GHEP mode** (the tool needs to instantiate providers
to compute responses). For flat-tree mode `-c` is **optional**: if given,
`applies_to_channels` patterns are read from it and used to skip booking
histograms for (dial × channel) combinations the dial doesn't apply to.

#### What gets plotted, what gets pruned

The keep/discard pipeline runs in this order:

1. **`applies_to_channels` filter** (fill time) — never book a histogram for
   (channel, dial) when the channel key doesn't match the dial's pattern list.
   Topology-level filter.
2. **`PruneEmpty`** — drop channels with zero events overall.
3. **Total aggregation** — emit `Total`, `Total_CC`, `Total_NC` only when at
   least one source channel had entries for the dial.
4. **`PruneTrivialParams`** — drop any (channel, dial) where every variation
   histogram matches CV bin-for-bin across every variable. Catches dials
   whose internal GENIE cuts (struck-nucleon PDG, FS pion multiplicity, W
   cuts, FS-hadron presence, etc.) reject every event in a channel even
   though the topology nominally matched.
5. **Print loop** — sub-pages of ≤6 variables per (channel, dial).

The `[INFO]` lines printed at runtime tell you what each stage did. On the
validation 3.06 ν<sub>μ</sub>+Ar40 1000-event sample the chain reduces from
1535 pages (no filter) → 411 pages (with all filters on).

---

### `DumpDialResponseNuSyst`

For each configured dial, samples `N` events per interaction channel and
plots `weight vs paramValue` curves. Uses
`GetEventResponse(evt, paramVals)` with `OverrideVariations` to scan an
arbitrary grid (use `--nsteps M` to resample the dial's declared variation
range with `M` evenly-spaced points).

```bash
DumpDialResponseNuSyst -c <config.fcl> -i <ghep.root> \
                       [-n <events-per-channel>] [--nsteps M] \
                       [-o <basename>] [-p MaCCQE,FormZone] \
                       [-N <NMax>] [--no-pdf | --no-root]
```

Useful when you want to see whether a dial moves an individual event monotonically,
whether it saturates, and where the curve flattens.

---

### `DeclaredDialTestNuSyst` (`nusyst inventory`)

Inventory dump — instantiates every provider declared in the config and
prints what each dial says about itself. No event input.

```bash
nusyst inventory                                      # auto-generate then inventory
nusyst inventory -c <config.fcl>                      # use existing fhicl
nusyst inventory --mode CCQE                          # filter to CCQE-bucket dials
nusyst inventory --mode FSI                           # filter to FSI dials
nusyst inventory -p MaCCQE,FormZone --verbose         # name-substring filter + detail
```

Run with no arguments: shells out to `GenerateAllDialsConfigNuSyst --mode all`,
writes the result to a temp file under `/tmp/nusyst_inventory_*.fcl`, runs
the inventory on it, and deletes the temp file on exit. Pass `--keep-tmp` to
preserve it for inspection. Pass `-c <config.fcl>` to inventory a config you
built yourself.

**`--mode <BUCKET>`** filters to dials whose provider maps to the requested
bucket. One bucket per provider — no "FSI is included in CCQE" overlap:

| `--mode` value | Providers shown |
|---|---|
| `CCQE` | `GENIEReWeight_CCQE`, `ZExpPCAWeighter_*` (b₁..b₄) |
| `CCRES` | `GENIEReWeight_CCRES` |
| `NCEL` | `GENIEReWeight_NCEL` |
| `NCRES` | `GENIEReWeight_NCRES` |
| `RES` | `GENIEReWeight_RES`, `ResIso_ResonanceIsolation` |
| `COH` | `GENIEReWeight_COH` |
| `DIS` | `GENIEReWeight_DIS` |
| `MEC` | `GENIEReWeight_MEC` |
| `SPP` | `GENIEReWeight_SPP` |
| `FSI` | `GENIEReWeight_FSI` |

The GENIE Reweight providers map by their `_<BUCKET>` suffix; non-RW
providers (`ZExpPCAWeighter`, `ResIso`) are hardcoded by tool_type in
`ModeForProvider()`.

`--verbose` prints the full `paramVariations`, validity range, and tool
options for each dial. `-p substr` filters dial names by substring (matches
multiple substrings if comma-separated). `-p` and `--mode` compose.

---

## Channel buckets and `applies_to_channels`

When `GenerateAllDialsConfigNuSyst` is invoked *without* `--single-instance`,
it groups every GENIE Reweight dial into one of these buckets:

| Bucket | Patterns | Dials |
|---|---|---|
| `CCQE` | `CC_*_*_QE` | `MaCCQE`, `CCQEPauliSupViaKF`, `CCQEMomDistroFGtoSF`, `RPA_CCQE`, `CoulombCCQE`, … |
| `CCRES` | `CC_*_*_RES` | `MaCCRES`, `MvCCRES` |
| `NCEL` | `NC_*_*_QE` | `MaNCEL`, `EtaNCEL` |
| `NCRES` | `NC_*_*_RES` | `MaNCRES`, `MvNCRES` |
| `RES` | `*_*_*_RES` | Delta internals (`RDecBR1eta`, `Theta_Delta2Npi`, …) — apply to CC+NC RES |
| `COH` | `*_*_*_COH` | `MaCOHpi`, `R0COHpi`, `NormCCCOH`, `NormNCCOH` |
| `DIS` | `*_*_*_DIS` | Bodek-Yang shape (`AhtBY`, `BhtBY`, `CV1uBY`, `CV2uBY`), AGKY hadronisation (`AGKYpT1pi`, `AGKYxF1pi`), `RnubarnuCC`, `NormDISCC` |
| `MEC` | `*_*_*_MEC` | `NormCCMEC`, `NormNCMEC`, `XSecShape_CCMEC`, `DecayAngMEC`, `FracDelta_CCMEC`, `FracPN_CCMEC` |
| `SPP` | `*_*_*_RES`, `*_*_*_DIS` | `NonRESBGv{p,n}{CC,NC}{1,2}pi`, `NonRESBGvbar...`, `MKSPP_ReWeight` |
| `FSI` | *(empty — applies to all)* | `FormZone`, `MFP_{pi,N}`, `Fr{CEx,Inel,Abs,PiProd,Elas}_{pi,N}` |

The channel key for an event is built by `nusyst::channel::MakeChannelKey`:

```
<CC|NC>_<nu-species>_<target-name>_<topology>
   e.g.  CC_numu_Ar40_QE
         NC_nuebar_C12_RES
```

`applies_to_channels` is a vector of fnmatch glob patterns (POSIX defaults,
so `*` matches across underscores, `?` matches a single character).
An empty list means "applies everywhere". The runtime filter lives in
`src/nusystematics/utility/channel_classification.hh`
(`nusyst::channel::MatchesAny`).

Downstream tools that honour `applies_to_channels`:

- `DumpConfiguredTweaksNuSyst` — skips evaluating non-matching providers
  per event, substituting a trivial response.
- `PlotSystVariationsNuSyst` — never books histograms for non-matching
  (dial × channel) combinations.

Pass `--single-instance` to disable the split — handy if you have a
downstream consumer that assumes one provider instance.

---

## CCQE axial-FF dials: raw vs. PCA-rotated

GENIE's CCQE axial form factor can be set up either as a 1-parameter dipole or
a multi-parameter Z-expansion. The Z-expansion fit (Meyer, Betancourt, Gran,
Hill — *PRD 93, 113015*) gives four coefficients `a₁..a₄` with a published
4×4 covariance matrix. There are two ways to vary them, and **a `--mode all`
config will emit both** so they can be compared:

### Raw `ZExpA1..4CCQE` + `ZNormCCQE` (provider: `GENIEReWeight_CCQE`)

Each dial shifts one coefficient independently:

```
a_i_new = a_i_default * (1 + twk_i * fracerr_zexp[i])
```

There is **no correlation matrix in this layer**. The per-coefficient
fractional errors come from GENIE's `CommonParam.xml`. Tossing all four at
±1σ ignores the (strong) anti-correlations among the `aₙ` and grossly
overcounts the constrained directions. Useful for one-at-a-time variation
studies; **wrong** for a covariance-aware envelope.

### PCA-rotated `b₁..b₄` (provider: `ZExpPCAWeighter`, fcl `zexpansion_weighter.ToolConfig.fcl`)

Diagonalises the published 4×4 covariance once with
`Eigen::SelfAdjointEigenSolver`:

```
Σ = V Λ Vᵀ                 // V = eigenvectors, Λ = diag(eigenvalues)
P = V · diag(√Λ)            // decorrelation transform (per-column √λₖ·vₖ)
```

The user-facing dials become `b₁..b₄` — **uncorrelated, unit-variance**
coordinates in the principal-component basis. At weight time:

```
a_shift     = P · b                                  // ChangeBasisBParams
a_for_genie = ((a_cv + a_shift)/a_genie_default - 1) / fracerr_zexp_genie
                                                     // ScaleAparamsforGenie
// → pushed into the underlying ZExpA1..4CCQE GReWeight dials.
```

So `b` is the basis you'd want to draw correlated toys from; the raw
`ZExpA*CCQE` are the basis the actual GENIE reweighter uses internally.

**For systematic envelopes / fits, use `b₁..b₄` only.** Both sets are
present in `--mode all` purely so you can compare per-event responses for
validation and to confirm the rotation is doing what it should. Enabling
both in a fit double-counts the same CCQE axial-FF physics in two bases.

The hardcoded covariance lives in `ZExpPCAWeighter_tool.cc`, in the
`namespace PRD_93_113015` block. The four `aₖ` central values
`(2.30, -0.60, -3.80, 2.30)` and per-coefficient errors
`(0.13, 1.0, 2.5, 2.7)` are from that paper's deuterium fit; the
`Covariance_Matrix` is the published `aᵢ`-basis covariance.
