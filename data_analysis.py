"""
data_analysis.py
================
Loads per-residue dihedral angle CSVs, computes histograms across simulation
repeats, calculates three distributional distance metrics, generates comparison
plots, and writes a ranked summary CSV.

This module is not intended to be run directly; call
``aggregate_and_plot_parallel()`` from ``dihedral_distribution_analysis.py``.

Output files
------------
png_output/residue_{id}_comparison.png
    Per-residue overlay plot: scatter of individual repeats + mean line.
txt_output/residue_{id}_output.txt
    Per-residue mean ± std histogram table plus metric values.
summary.csv
    All residues ranked by TVD (descending), with TVD, JSD, and Wasserstein.
output.txt
    Master processing log with timestamps.

Symmetric folding
-----------------
When ``fold_symmetric=True`` is passed, raw angles in [−180°, +180°] are
transformed to [0°, 180°) via ``angle % 180`` before histogramming.  This is
appropriate for residues whose terminal group has C2 symmetry (e.g. PHE/TYR
chi2, ASP/GLU chi2), where the two halves of the full dihedral range are
physically equivalent.  The raw CSV files are never modified.
"""

import os
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool

from scipy.spatial.distance import cityblock, jensenshannon
from scipy.stats import wasserstein_distance

# ---------------------------------------------------------------------------
# Output directories (created at import time so workers can always write)
# ---------------------------------------------------------------------------
os.makedirs("png_output", exist_ok=True)
os.makedirs("txt_output", exist_ok=True)
os.makedirs("csv_output", exist_ok=True)


# ---------------------------------------------------------------------------
# Per-residue worker
# ---------------------------------------------------------------------------

def process_residue(args: tuple) -> tuple:
    """
    Compute histograms, distance metrics, and generate a comparison plot for
    one residue.

    Parameters
    ----------
    args : tuple
        (resid_key_and_files, degrees_range, simulation_1_label,
         simulation_2_label, dihedral_label, fold_symmetric)
        where ``resid_key_and_files`` is (resid_str, file_dict).

    Returns
    -------
    (log_message, result_tuple | None)
        result_tuple is (resid, tvd, jsd, wd) on success, else None.
    """
    resid_key_and_files, degrees_range, simulation_1_label, simulation_2_label, dihedral_label, fold_symmetric = args
    resid, file_dict = resid_key_and_files

    # Bin edges depend on whether we are folding to [0°, 180°) or keeping [−180°, +180°]
    if fold_symmetric:
        bin_edges = np.arange(0, 180 + degrees_range, degrees_range)
    else:
        bin_edges = np.arange(-180, 180 + degrees_range, degrees_range)

    sim1_files = file_dict[resid].get("1", [])
    sim2_files = file_dict[resid].get("2", [])

    if not sim1_files or not sim2_files:
        return f"[SKIP] {resid} – missing files for one or both simulations", None

    # ------------------------------------------------------------------
    def load_histograms(file_list: list):
        """Return an (n_repeats × n_bins) array of raw counts, or a str on error."""
        histograms = []
        for f in file_list:
            try:
                data = pd.read_csv(f, header=None).iloc[:, 0].dropna()
                if fold_symmetric:
                    # Fold [−180°, +180°] → [0°, 180°) using C2 symmetry:
                    #   angle % 180  maps −170° → 10°, −90° → 90°, 170° → 170°
                    # This collapses the two physically equivalent halves of
                    # symmetric groups (e.g. PHE/TYR chi2, ASP/GLU chi2) onto
                    # a single [0°, 180°] axis.
                    data = data % 180
                hist, _ = np.histogram(data, bins=bin_edges)
                histograms.append(hist)
            except Exception as exc:
                return f"[ERROR] Could not load {f}: {exc}"
        return np.array(histograms)
    # ------------------------------------------------------------------

    sim1_hists = load_histograms(sim1_files)
    sim2_hists = load_histograms(sim2_files)

    if isinstance(sim1_hists, str):
        return sim1_hists, None
    if isinstance(sim2_hists, str):
        return sim2_hists, None
    if sim1_hists.size == 0 or sim2_hists.size == 0:
        return f"[SKIP] {resid} – empty histogram array", None

    # Mean and standard deviation across repeats
    mean1, std1 = sim1_hists.mean(axis=0), sim1_hists.std(axis=0)
    mean2, std2 = sim2_hists.mean(axis=0), sim2_hists.std(axis=0)

    if mean1.sum() == 0 or mean2.sum() == 0:
        return f"[SKIP] {resid} – zero-sum histogram (no data in range)", None

    # Normalise to probability distributions
    p = mean1 / mean1.sum()
    q = mean2 / mean2.sum()

    # Distance metrics
    tvd = 0.5 * cityblock(p, q)       # Total Variation Distance       → [0, 1]
    jsd = jensenshannon(p, q)          # Jensen–Shannon Divergence (√)  → [0, 1]
    wd  = wasserstein_distance(p, q)   # Wasserstein / Earth Mover's    → [0, ∞)

    # ------------------------------------------------------------------
    # Text report
    # ------------------------------------------------------------------
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    fold_note = " [symmetric fold: 0–180°]" if fold_symmetric else ""
    lines = [
        f"Residue {resid}{fold_note}",
        f"  TVD:              {tvd:.4f}",
        f"  Jensen-Shannon:   {jsd:.4f}",
        f"  Wasserstein:      {wd:.4f}",
        f"  Repeats (sim 1):  {len(sim1_files)}",
        f"  Repeats (sim 2):  {len(sim2_files)}",
        "",
        f"  {simulation_1_label} – mean ± std per bin",
    ]
    for i, (m, s) in enumerate(zip(mean1, std1)):
        lines.append(f"    {bin_edges[i]:.0f}° to {bin_edges[i+1]:.0f}°:  {m:.2f} ± {s:.2f}")
    lines.append(f"\n  {simulation_2_label} – mean ± std per bin")
    for i, (m, s) in enumerate(zip(mean2, std2)):
        lines.append(f"    {bin_edges[i]:.0f}° to {bin_edges[i+1]:.0f}°:  {m:.2f} ± {s:.2f}")
    lines.append("")

    txt_path = os.path.join("txt_output", f"residue_{resid}_output.txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines))

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))

    # Individual repeats as scatter (show variance)
    for hist in sim1_hists:
        ax.scatter(bin_centers, hist, alpha=0.35, color="C0", s=25)
    for hist in sim2_hists:
        ax.scatter(bin_centers, hist, alpha=0.35, color="C1", s=25)

    # Mean lines
    ax.plot(bin_centers, mean1, "-o", color="C0", label=simulation_1_label, linewidth=1.8)
    ax.plot(bin_centers, mean2, "-o", color="C1", label=simulation_2_label, linewidth=1.8)

    if fold_symmetric:
        ax.set_xticks(np.arange(0, 181, 30))
        ax.set_xlabel(f"{dihedral_label} dihedral (°)  [folded: 0–180°]", fontsize=22)
    else:
        ax.set_xticks(np.arange(-180, 181, 30))
        ax.set_xlabel(f"{dihedral_label} dihedral (°)", fontsize=22)

    ax.tick_params(labelsize=18)
    ax.set_ylabel("Frequency", fontsize=22)
    ax.set_title(
        f"{dihedral_label} – {resid}{fold_note}   "
        f"[TVD={tvd:.3f}  JSD={jsd:.3f}  WD={wd:.3f}]",
        fontsize=16,
    )
    ax.legend(fontsize=14)
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join("png_output", f"residue_{resid}_comparison.png"), dpi=150)
    plt.close(fig)

    return f"[DONE] {resid}", (resid, tvd, jsd, wd)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def aggregate_and_plot_parallel(
    degrees_range: int,
    simulation_1_label: str,
    simulation_2_label: str,
    dihedral_label: str,
    fold_symmetric: bool = False,
) -> None:
    """
    Discover per-residue CSV files, dispatch per-residue analysis in parallel,
    and write summary outputs.

    Parameters
    ----------
    degrees_range : int
        Histogram bin width in degrees.
    simulation_1_label : str
        Human-readable label for simulation 1 (used in plots and CSV header).
    simulation_2_label : str
        Human-readable label for simulation 2.
    dihedral_label : str
        "chi1" or "chi2"; used for file globbing and axis labels.
    fold_symmetric : bool, optional
        If True, angles are folded from [−180°, +180°] to [0°, 180°) before
        histogramming, by applying ``angle % 180``.  Appropriate for residues
        whose terminal group has C2 symmetry (PHE/TYR chi2, ASP/GLU chi2).
        Default is False.
    """
    # Discover and organise CSV files
    file_dict: dict = {}
    pattern = f"simulation_*_repeat_*_{dihedral_label}_resid_*.csv"
    for csv_file in glob.glob(pattern):
        parts = csv_file.split("_")
        if len(parts) < 7:
            print(f"[WARN] Unexpected filename format, skipping: {csv_file}")
            continue
        sim_number = parts[1]                    # "1" or "2"
        residue_id = parts[6].split(".")[0]      # e.g. "VAL62"
        file_dict.setdefault(residue_id, {"1": [], "2": []})
        if sim_number in ("1", "2"):
            file_dict[residue_id][sim_number].append(csv_file)

    if not file_dict:
        print(
            f"[ERROR] No CSV files matched pattern '{pattern}'.\n"
            "        Check that extraction has been run and that DIHEDRAL_LABEL is correct."
        )
        return

    print(f"[INFO] Found {len(file_dict)} residue(s) to analyse.")

    # Build task list and run in parallel
    task_args = [
        ((resid, file_dict), degrees_range, simulation_1_label, simulation_2_label, dihedral_label, fold_symmetric)
        for resid in file_dict
    ]

    with Pool() as pool:
        results = pool.map(process_residue, task_args)

    # ------------------------------------------------------------------
    # Master log
    # ------------------------------------------------------------------
    fold_note = " (symmetric fold: 0–180°)" if fold_symmetric else ""
    with open("output.txt", "w") as fh:
        fh.write(f"Dihedral distribution analysis  –  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"Dihedral      : {dihedral_label}{fold_note}\n")
        fh.write(f"Sim 1         : {simulation_1_label}\n")
        fh.write(f"Sim 2         : {simulation_2_label}\n\n")
        for log_msg, _ in results:
            fh.write(log_msg + "\n")

    # ------------------------------------------------------------------
    # Summary CSV
    # ------------------------------------------------------------------
    summary_rows = []
    for _, result in results:
        if result is not None:
            resid, tvd, jsd, wd = result
            summary_rows.append({
                "residue":         resid,
                "TVD":             round(tvd, 4),
                "Jensen-Shannon":  round(jsd, 4),
                "Wasserstein":     round(wd,  4),
            })

    if not summary_rows:
        print("[WARNING] No residues were successfully processed. Check output.txt for details.")
        return

    summary_df = (
        pd.DataFrame(summary_rows)
        .sort_values(by="TVD", ascending=False)
        .reset_index(drop=True)
    )
    summary_df.to_csv("summary.csv", index=False)

    # ------------------------------------------------------------------
    # Global statistics
    # ------------------------------------------------------------------
    done_count = sum(1 for m, _ in results if m.startswith("[DONE]"))
    skip_count = sum(1 for m, _ in results if m.startswith("[SKIP]"))

    print(f"\n{'='*50}")
    print(f"  ANALYSIS COMPLETE")
    print(f"  Residues processed : {done_count}")
    print(f"  Residues skipped   : {skip_count}")
    print(f"{'='*50}")
    print(f"  Mean TVD             : {summary_df['TVD'].mean():.4f}")
    print(f"  Mean Jensen-Shannon  : {summary_df['Jensen-Shannon'].mean():.4f}")
    print(f"  Mean Wasserstein     : {summary_df['Wasserstein'].mean():.4f}")
    print(f"{'='*50}")
    print(f"\n  Top 10 residues by TVD:\n")
    print(summary_df.head(10).to_string(index=False))
    print(f"\n  Full results → summary.csv")
