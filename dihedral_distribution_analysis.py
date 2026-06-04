"""
dihedral_distribution_analysis.py
==================================
Entry point and user configuration for the dihedral distribution pipeline.

Edit the USER SETTINGS section below, then run:
    python dihedral_distribution_analysis.py

Two modes are available via the RUN_* flags:
    EXTRACT  – parse MD trajectories and write per-residue angle CSVs
    ANALYSE  – load those CSVs, compute statistics, and generate plots
Both are enabled by default; disable either if you have already run that step.
"""

import glob
import shutil

from dihedral_definitions import DihedralDefinitions_chi1, DihedralDefinitions_chi2
from data_extraction import run_parallel_extraction
from data_analysis import aggregate_and_plot_parallel

# =============================================================================
# USER SETTINGS  – edit this section only
# =============================================================================

DIHEDRAL_LABEL = "chi1"
DEGREES_RANGE = 10

SIMULATION_1_LABEL = "BNB"
SIMULATION_2_LABEL = "BNO"

SEGID_OF_INTEREST = "P0"   # IMPORTANT: If you prefer chainid change "segid" to "chainid" in all .py files!

FOLD_SYMMETRIC = False

# --- Pipeline stages to run ---
RUN_EXTRACTION = True   # Set False to skip trajectory parsing (re-use existing CSVs)
RUN_ANALYSIS   = True   # Set False to skip statistical analysis / plotting

# --- Simulation file paths --- Add as many as you like
# Format: { "simulation_1": [(traj1, top1), (traj2, top2), ...],
#           "simulation_2": [(traj1, top1), ...] }
# Keys MUST be "simulation_1" and "simulation_2".
SIMULATIONS = {
    "simulation_1": [
        (
            "system_x/repeat_1/trajectory_repeat_1.xtc",
            "system_x/repeat_1/structure.pdb",
        ),
        (
            "system_x/repeat_2/trajectory_repeat_3.xtc",
            "system_x/repeat_2/structure.pdb",
        ),
            "system_x/repeat_3/trajectory_repeat_3.xtc",
            "system_x/repeat_3/structure.pdb",
        ),

    ],
    "simulation_2": [
        (
            "system_y/repeat_1/trajectory_repeat_1.xtc",
            "system_y/repeat_1/structure.pdb",
        ),
        (
            "system_y/repeat_2/trajectory_repeat_3.xtc",
            "system_y/repeat_2/structure.pdb",
        ),
            "system_y/repeat_3/trajectory_repeat_3.xtc",
            "system_y/repeat_3/structure.pdb",
        ),
    ],
}

# =============================================================================
# END OF USER SETTINGS
# =============================================================================

# Select the correct dihedral definitions based on user choice
if DIHEDRAL_LABEL == "chi1":
    dihedral_definitions = DihedralDefinitions_chi1
elif DIHEDRAL_LABEL == "chi2":
    dihedral_definitions = DihedralDefinitions_chi2
else:
    raise ValueError(f"Unknown DIHEDRAL_LABEL '{DIHEDRAL_LABEL}'. Choose 'chi1' or 'chi2'.")


if __name__ == "__main__":

    if RUN_EXTRACTION:
        print(f"\n[INFO] Starting dihedral extraction ({DIHEDRAL_LABEL}) ...")
        run_parallel_extraction(
            simulations_dict=SIMULATIONS,
            segid=SEGID_OF_INTEREST,
            dihedral_label=DIHEDRAL_LABEL,
            dihedral_definitions=dihedral_definitions,
        )

    if RUN_ANALYSIS:
        fold_note = " [symmetric fold: 0–180°]" if FOLD_SYMMETRIC else ""
        print(f"\n[INFO] Starting distribution analysis ({DIHEDRAL_LABEL}{fold_note}) ...")
        aggregate_and_plot_parallel(
            degrees_range=int(DEGREES_RANGE),
            simulation_1_label=SIMULATION_1_LABEL,
            simulation_2_label=SIMULATION_2_LABEL,
            dihedral_label=DIHEDRAL_LABEL,
            fold_symmetric=FOLD_SYMMETRIC,
        )

    # Archive raw angle CSVs (everything except summary.csv)
    archived = 0
    for csv_file in glob.glob("*.csv"):
        if csv_file != "summary.csv":
            shutil.move(csv_file, "csv_output/")
            archived += 1
    if archived:
        print(f"\n[INFO] Moved {archived} raw CSV file(s) to csv_output/")
