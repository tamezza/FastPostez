# fastpostez — Post-processing of EasyJet Ntuples for b-JES In-Situ Calibration

## Overview

Super-fast post-processing framework for [easyjet](https://gitlab.cern.ch/easyjet/easyjet) ntuples. It applies physics selections, computes derived quantities, and produces output ntuples used for b-JES and b-JER in-situ calibration studies for $Z \rightarrow bb + \gamma$ and $Z \rightarrow bb + \mathrm{jets}$ topologies.

---

## Setup

```bash
source setup.sh     # Set up ROOT and compiler environment
./compile.sh        # Configure and build with CMake
```

`setup.sh` sets up the ATLAS / ROOT environment used on `lpnatlas02` cluster; `compile.sh` runs CMake and builds everything into `build/`. After a successful build, the main executable `run_bjes` is on your `PATH`.

---

## Workflow

The post-processing for each channel is two steps:

1. **`run_bjes`** — read the easyjet input, apply selection + derive quantities, and write per-sample skimmed ntuples.
2. **`merge_zbbj*.sh` / `merge_zbby*.sh`** — hadd the per-sample outputs into the final per-process ntuples used by the analysis.

General form:

```bash
run_bjes -c <config.json> -i <input_dir> -o <output_dir>
./merge_<channel>[_run3].sh <run_bjes_output_dir> <final_merged_dir>
```

The configs (`zbbj_bjes_run{2,3}_config.json`, `zbby_bjes_run{2,3}_config.json`) define the selection for each channel and data-taking period.

---

## Directory Structure

```
fastpostez/
├── bjes/                        # Core analysis code
│   ├── bjes.cpp                 # Main analysis logic (RDataFrame)
│   ├── bjes.h
│   └── CMakeLists.txt
├── core/                        # Shared utilities
│   ├── input_handler.cpp/h      # File discovery and TChain building
│   ├── weights.cpp/h            # Event weight computation
│   ├── output_manager.cpp/h     # Output directory management
│   ├── config_bjes.h            # Config struct definition
│   └── gn2x.cpp/h               # GN2X b-tagging WP handler
├── utils/
│   ├── get_metadata.py          # Fetch xsec/filtEff from AMI
│   └── setup.sh                 # ATLAS environment setup
├── data/
│   ├── metadata.txt             # xsec, genFiltEff, kFactor per DSID
│   └── Xbb_lookup_table_*.json  # GN2X flat-mass WP lookup table
├── zbbj_bjes_run3_config.json   # Config for Z(bb)+jets Run 3
├── zbby_bjes_run3_config.json   # Config for Z(bb)+γ Run 3
├── zbbj_bjes_run2_config.json   # Config for Z(bb)+jets Run 2
├── zbby_bjes_run2_config.json   # Config for Z(bb)+γ Run 2
├── compile.sh                   # Build script
└── setup.sh                     # Environment setup
```

---

## Configuration

Each analysis requires a JSON config file. Example (`zbby_bjes_run3_config.json`):

```json
{
    "ntuple": {
        "tree_name": "AnalysisMiniTree"
    },
    "analysis": {
        "dsids": [545023, 545024, 701286, 701287, 701288],
        "years": [22, 23, 24],
        "trigger_map": {
            "2022": ["HLT_g140_loose_L1EM22VHI"],
            "2023": ["HLT_g140_loose_L1EM22VHI", "HLT_g140_loose_L1eEM26M"],
            "2024": ["HLT_g140_loose_L1eEM26M",  "HLT_g140_loose_L1EM22VHI"]
        },
        "jet_jvt_branch": "recojet_antikt4PFlow_jvt_selection_NOSYS",
        "metadata": "data/metadata.txt",
        "flatmass": "data/Xbb_lookup_table_prelim_Oct30_2024.json",
        "wps": ["FlatMassQCDEff_0p46", "FlatMassQCDEff_0p58"],
        "analysis": "zbby"
    }
}
```

### Key config fields

| Field | Description |
|---|---|
| `tree_name` | Input tree name (typically `AnalysisMiniTree`) |
| `dsids` | List of MC DSIDs to process |
| `years` | Data-taking years for data samples (22, 23, 24) |
| `trigger_map` | Per-year trigger requirements |
| `jet_jvt_branch` | JVT branch name — differs between samples (see below) |
| `metadata` | Path to xsec/filtEff/kFactor file |
| `flatmass` | Path to GN2X flat-mass WP lookup table |
| `wps` | List of GN2X working points to evaluate |
| `analysis` | Analysis type: `"zbbj"` or `"zbby"` |

### JVT branch name per sample type

| Sample | JVT branch |
|---|---|
| zbby (Run 3) | `recojet_antikt4PFlow_jvt_selection_NOSYS` |
| zbbj (Run 3) | `recojet_antikt4PFlow_Jvt_NOSYS` |

The code resolves this automatically via `get_branch_name()` - the config field is optional but recommended for clarity.

---

## Running

```bash
run_bjes -c <config.json> -i <input_folder> -o <output_folder>
```

### Options

| Flag | Description |
|---|---|
| `-c` | Config JSON file |
| `-i` | Input folder containing per-DSID subdirectories |
| `-o` | Output folder (created automatically) |
| `-n` | Number of threads for RDataFrame (default: 1) |

### Examples

```bash
# Z(bb)+gamma Run 3, single-threaded
run_bjes -c zbby_bjes_run3_config.json \
         -i /data/atlas/tamezza/BJetCalib_Run3/easyjet_ntuples/samples_zbby/ybjes_run3 \
         -o /data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/JetCalibxbb_Zbby_GN2x_run3

# Z(bb)+jets Run 3, 8 threads
run_bjes -c zbbj_bjes_run3_config.json \
         -i /data/atlas/tamezza/BJetCalib_Run3/easyjet_ntuples/samples_zbbj/bjes_run3 \
         -o /data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/JetCalibxbb_Zbbj_GN2x_run3 \
         -n 8
```

---

## Input File Structure

`fastpostez` expects input files organized by DSID. For each DSID, it searches the input folder for subdirectories matching the DSID number and collects all `*.output-tree.root` files. If a `merged.root` exists it is used directly; otherwise individual files are chained.

**Important:** If `merged.root` exists but is stale (produced before new branches were added), delete it and rebuild:

```bash
# Delete all stale merged files
find <input_folder> -name "merged.root" -delete

# Rebuild merged files for all DSIDs
for dir in <input_folder>/*/; do
    hadd -f "${dir}merged.root" "${dir}"user.tamezza.*.output-tree.root
done
```

---

## Output Structure

```
<output_folder>/
├── ntuples/
│   ├── ntuples_545023.root    # Per-DSID output ntuple (tree: "nominal")
│   ├── ntuples_545024.root
│   ├── ntuples_22.root        # Per-year data ntuple
│   ├── ntuples_23.root
│   └── ntuples_24.root
└── histograms/
    └── hists_<sample>.root    # Cutflow and summary histograms
```

### Output branches (common to all samples)

| Branch | Type | Description |
|---|---|---|
| `ljet_pt/eta/phi/m/e` | `RVecF` | Large-R jet kinematics (GeV) |
| `jet_pt/eta/phi/m/e` | `RVecF` | Small-R jet kinematics (GeV) |
| `jet_jvt` | `RVecF` | Small-R jet JVT score |
| `ljet_bJR10v00_pt/mass` | `RVecF` | bJR10v00 variant kinematics |
| `ljet_bJR10v01_pt/mass` | `RVecF` | bJR10v01 variant kinematics |
| `dhbb` | `RVecF` | GN2X discriminant score |
| `pass_FlatMassQCDEff_*` | `RVecB` | GN2X WP pass/fail flags |
| `total_weight` | `double` | Full event weight (xsec × lumi × genFiltEff × pileup) |
| `n_ljets`, `n_jets` | `int` | Jet multiplicities |

### MC-only branches

| Branch | Type | Description |
|---|---|---|
| `ljet_truth_label` | `RVecI` | Large-R jet truth label |
| `jet_parton_truth_label` | `RVecI` | Small-R jet parton label |
| `jet_eff_jvt` | `float` | JVT efficiency scale factor |
| `ljet_n_bhadrons` | `RVecI` | Ghost-associated b-hadron count |
| `ljet_n_chadrons` | `RVecI` | Ghost-associated c-hadron count |
| `matched_tjet_pt/m/eta/phi` | `RVecF` | Truth-matched large-R jet kinematics (GeV) |
| `matched_tjet_dR` | `RVecF` | ΔR between reco and matched truth jet |

### zbby-only branches

| Branch | Type | Description |
|---|---|---|
| `ph_pt/eta/phi/e` | `RVecF` | Photon kinematics (GeV) |
| `met_met`, `met_phi` | `float` | Missing ET |

---

## Metadata

The metadata file provides cross-sections, generator filter efficiencies, and k-factors for event weight normalization. Format:

```
dataset_number/I:crossSection_pb/D:genFiltEff/D:kFactor/D:generator_name/C
545023  <xsec>  <filtEff>  1.0  PowhegPy8
701286  <xsec>  <filtEff>  1.0  Sherpa(v.2.2.14)
```

### Fetching metadata from AMI

First set up the ATLAS environment:

```bash
source setup.sh
# or
lsetup "asetup AthAnalysis,25.2.55" pyAMI
```

Then run:

```bash
cd utils/
python get_metadata.py --campaign mc23 --output-mc23 ../data/metadata_zbby_run3.txt
```

**Note:** Different sample types use different reconstruction tags (`rtag`). The zbby signal samples (545023, 545024, 701286–701288) use `r15530`, while other Run 3 samples use `r15224_r15225`. The script handles this by calling `get_metadata()` separately for each group.

Append new entries to the main metadata file:

```bash
tail -n +2 ../data/metadata_zbby_run3.txt >> ../data/metadata.txt
tail -n +2 ../data/metadata_zbby_signal.txt >> ../data/metadata.txt
```

---

## Known Branch Name Differences

Some branches have different names across sample types. The code resolves these automatically using `get_branch_name()` with a fallback list:

| Logical name | zbby branch | zbbj branch |
|---|---|---|
| `jet_jvt` | `recojet_antikt4PFlow_jvt_selection_NOSYS` | `recojet_antikt4PFlow_Jvt_NOSYS` |
| `jet_eff_jvt` | `jvt_effSF_NOSYS` | `recojet_antikt4PFlow_jvt_effSF_NOSYS` |
| `ljet_bJR10v00_mass` | `recojet_antikt10UFO_bJR10v00Ext_mass_NOSYS` | `recojet_antikt10UFO_bJR10v00_mass_NOSYS` |

---

## Common Issues

**`Failed to create valid TChain from input files`**
All files for a DSID have 0 entries. Either the grid jobs failed silently or the files are corrupt. Check with:
```bash
for f in <dsid_dir>/*.root; do
    root -l -b -q -e "TFile f(\"$f\"); auto t=(TTree*)f.Get(\"AnalysisMiniTree\"); std::cout<<(t?t->GetEntries():0)<<\"\n\";" 2>/dev/null
done
```

**`use of undeclared identifier '<branch_name>'`**
A branch required by the config does not exist in the input file. Usually caused by a stale `merged.root`. Delete it and rebuild (see above).

**Negative event weights**
The DSID is missing from `metadata.txt`. Add it using `get_metadata.py` and rerun.

**`getMetadata.py: command not found`**
The ATLAS environment is not set up. Run `source setup.sh` or `lsetup "asetup AthAnalysis,25.2.55" pyAMI` before running the script.

---
## Authors

- Thiziri Amezza (thiziri.amezza@cern.ch)
