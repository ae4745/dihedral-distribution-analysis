# Dihedral Distribution Analysis

A Python pipeline for comparing side-chain dihedral angle distributions (χ1 / χ2) across molecular dynamics simulation ensembles.

Given two sets of simulation repeats, the pipeline:
1. Extracts per-residue dihedral angles from MD trajectories (via MDAnalysis)
2. Computes histograms across repeats and compares them using three complementary statistical distance metrics
3. Generates per-residue comparison plots and a ranked summary CSV


# Background

### Side-chain dihedral angles

χ1 and χ2 are dihedral angles that describe side-chain rotameric states:

| Angle | Atoms | Measures |
|-------|-------|----------|
| χ1 | N – CA – CB – Xγ | Rotation about the CA–CB bond |
| χ2 | CA – CB – CG – Xδ | Rotation about the CB–CG bond |

Comparing their distributions between two simulation conditions (e.g. wild-type vs mutant, ligand-bound vs apo) reveals conformational changes at single-residue resolution, without reducing the data to a single rotamer label.

### Why use distributions and not just means?

Side-chain dihedrals are angular and multimodal (rotamer wells at roughly −60°, 60°, 180°). Comparing means is therefore misleading. This pipeline compares the full probability distributions using three metrics (TVD, Jensen-Shannon, Wasserstein) that are sensitive to different aspects of distributional change.

---

## Repository Structure

```
.
├── dihedral_distribution_analysis.py   ← START HERE: user settings & entry point
├── data_extraction.py                  # MD trajectory parsing (MDAnalysis)
├── data_analysis.py                    # Histogram comparison & statistics
├── dihedral_definitions.py             # χ1 / χ2 atom selections per residue type
├── requirements.txt
├── README.md
│
├── png_output/                         # Per-residue comparison plots   (auto-created)
├── txt_output/                         # Per-residue text reports        (auto-created)
└── csv_output/                         # Raw per-frame angle CSVs        (auto-created)
```

**You should only need to edit `dihedral_distribution_analysis.py`.**

---

## Dependencies

| Package | Minimum version |
|---------|----------------|
| MDAnalysis | 2.0 |
| NumPy | 1.21 |
| pandas | 1.3 |
| matplotlib | 3.4 | 
| SciPy | 1.7 | 
| tqdm | 4.62 |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ae4745/dihedral-distribution-analysis.git
cd dihedral-distribution-analysis
```

### 2. Create and activate a conda environment (recommended)

```bash
conda create -n dihedrals
conda activate dihedrals
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```


---

## Quick Start

1. Open `dihedral_distribution_analysis.py` and edit the **USER SETTINGS** section:
   - Set `DIHEDRAL_LABEL` to `"chi1"` or `"chi2"`
   - Set `SIMULATION_1_LABEL` and `SIMULATION_2_LABEL` to descriptive names
   - Set `SEGID_OF_INTEREST` to your chain identifier (e.g. `"P0"`, `"A"`)
      **Important! if you prefer a chain ID instead, you should replace "segid" with "chainid" in all .py files**
   - Populate `SIMULATIONS` with your trajectory and topology paths

2. Run the pipeline:

```bash
python dihedral_distribution_analysis.py
```

3. Check the outputs:
   - `summary.csv` — all residues ranked by TVD
   - `png_output/` — per-residue comparison plots
   - `txt_output/` — per-residue detailed reports

---

## Configuration Guide

All settings live in the `USER SETTINGS` block at the top of `dihedral_distribution_analysis.py`.

### `DIHEDRAL_LABEL`
```python
DIHEDRAL_LABEL = "chi1"   # or "chi2"
```
Selects which dihedral angle to analyse. This controls which atom selections are used and which CSV files are matched during analysis.

### `FOLD_SYMMETRIC`
```python
FOLD_SYMMETRIC = False   # or True
```
When `True`, raw dihedral angles in [−180°, +180°] are transformed to [0°, 180°) before histogramming, using:

```
angle_folded = angle % 180
```

This is the correct treatment for any terminal group with **C2 symmetry**, where a rotation of X° and X+180° are physically indistinguishable:

| Residue | Dihedral | Reason |
|---------|----------|--------|
| PHE, TYR | χ2 | Phenyl / tyrosyl ring; the two faces are equivalent |
| ASP, GLU | χ2 | Carboxylate; the two oxygens are equivalent |
| TRP | χ2 | Indole ring has approximate C2 symmetry about CB–CG |

Without folding, a rotamer that sits at −80° in one simulation and +100° in another will appear as two distinct populations even though they represent the same physical state. Folding collapses both onto 100° (since −80 % 180 = 100), and the histogram — and all three distance metrics — will correctly show no difference.

**When *not* to fold:** χ1 for any residue, and χ2 for residues without symmetric end-groups (ARG, GLN, LYS, MET, etc.). For mixed datasets, run the analysis twice (once with, once without) and compare.

> **Note:** Folding is applied only at the analysis stage. The raw CSVs always retain the original [−180°, +180°] values, so you can re-run with a different setting without re-extracting.

### `DEGREES_RANGE`
```python
DEGREES_RANGE = 10   # degrees per histogram bin
```
Bin width for the dihedral histograms. Smaller values give finer resolution but noisier distributions; 5–15° is typical for MD data.

### `SIMULATION_1_LABEL` / `SIMULATION_2_LABEL`
```python
SIMULATION_1_LABEL = "BNB"
SIMULATION_2_LABEL = "BNO"
```
Human-readable labels for the two conditions. Used in plot legends and the summary CSV header.

### `SEGID_OF_INTEREST`
```python
SEGID_OF_INTEREST = "P0"
```
The segment / chain identifier for the protein of interest in your topology.

### `RUN_EXTRACTION` / `RUN_ANALYSIS`
```python
RUN_EXTRACTION = True   # parse trajectories and write angle CSVs
RUN_ANALYSIS   = True   # load CSVs, compute metrics, generate plots
```
Set `RUN_EXTRACTION = False` if you have already extracted angles and only want to re-run the analysis (e.g. with a different bin width or labels). The raw CSVs are archived to `csv_output/` after each run; move them back to the working directory first.

### `SIMULATIONS`
```python
SIMULATIONS = {
    "simulation_1": [
        ("path/to/traj_rep1.xtc", "path/to/topology.pdb"),
        ("path/to/traj_rep2.xtc", "path/to/topology.pdb"),
        ...
    ],
    "simulation_2": [
        ...
    ],
}
```
Keys **must** be `"simulation_1"` and `"simulation_2"`. Each value is a list of `(trajectory, topology)` tuples, one per simulation repeat. The repeat index is assigned by list position (1-based). Any number of repeats is supported; unequal repeat counts between conditions are handled gracefully.

---

## Understanding the Outputs

### `summary.csv`

The primary output. One row per residue, sorted by TVD descending:

| Column | Description |
|--------|-------------|
| `residue` | Residue name + number (e.g. `VAL62`) |
| `TVD` | Total Variation Distance (0–1) |
| `Jensen-Shannon` | Jensen–Shannon Divergence (0–1) |
| `Wasserstein` | Wasserstein / Earth Mover's Distance |

### `png_output/residue_{id}_comparison.png`

One plot per residue showing:
- **Scatter points** — histogram counts for each individual repeat (shows inter-repeat variance)
- **Lines with markers** — mean histogram across repeats for each condition
- **Title** — includes all three metric values for quick visual triage

### `txt_output/residue_{id}_output.txt`

Per-residue text report with metric values and the full mean ± std histogram table for each condition.

### `csv_output/`

Raw per-frame dihedral angle values for every residue and repeat. Archived here after analysis. Filename format:
```
simulation_{N}_repeat_{R}_{dihedral}_{resname}{resid}.csv
```

### `output.txt`

Master processing log with a status line (`[DONE]`, `[SKIP]`, or `[FAIL]`) for every residue in every repeat.

---

## Interpreting the Metrics

### Total Variation Distance (TVD)

> TVD = ½ · Σ |p(x) − q(x)|

The maximum probability mass that is in a different bin between the two distributions. Ranges from 0 (identical) to 1 (no overlap).

**Practical thresholds (approximate):**

| TVD | Interpretation |
|-----|---------------|
| < 0.10 | Distributions are essentially the same |
| 0.10 – 0.20 | Subtle shift; may or may not be biologically meaningful |
| 0.20 – 0.35 | Moderate shift; worth inspecting the plot |
| > 0.35 | Large conformational change; high-priority residue |

### Jensen–Shannon Divergence (JSD)

A symmetrised and smoothed version of KL divergence. More robust than TVD when distributions have zero-probability bins (e.g. a rotamer that is completely absent in one condition). The square-root form (returned by SciPy) ranges from 0 to 1.

### Wasserstein Distance (Earth Mover's Distance)

Accounts for the geometry of the dihedral axis: a distribution shifted by 30° scores differently from one where probability mass teleports to the opposite side of the plot. Useful for detecting rotamer-well shifts vs. rotamer-population changes. Has no fixed upper bound; interpret relative to other residues in the same dataset.

### Using multiple metrics together

No single metric tells the whole story. A recommended workflow:

1. **Sort by TVD** — use as the primary ranking
2. **Cross-check with JSD** — if TVD is high but JSD is low, the shift may be concentrated in bins where one condition has very few counts (inspect the plot)
3. **Use Wasserstein** to distinguish *shifts* (the peak moved) from *redistributions* (mass moved between non-adjacent wells)

---

## Extending the Code

### Adding a new residue type

Edit `dihedral_definitions.py`:

```python
class DihedralDefinitions_chi1:
    Resids = {
        ...
        "MSE": "name N CA CB CG",   # selenomethionine example
    }
```

Use `None` if the dihedral does not exist for that residue.

### Analysing more than two conditions

The current design compares exactly two conditions. To extend:
- Add `"simulation_3"` etc. to `SIMULATIONS`
- Update `data_extraction.py` to handle extra simulation numbers
- Update `data_analysis.py` to load and compare additional histogram sets

