"""
data_extraction.py
==================
Extracts chi1 or chi2 dihedral angles from MD trajectories using MDAnalysis.

For each simulation repeat, each eligible residue is processed in parallel,
writing one CSV of angle values (one per frame) to the working directory.

Output filename format:
    simulation_{N}_repeat_{R}_{dihedral}_{resname}{resid}.csv

This module is not intended to be run directly; import and call
``run_parallel_extraction()`` from ``dihedral_distribution_analysis.py``.
"""

import os

import numpy as np
import MDAnalysis as mda
from multiprocessing import Pool, cpu_count
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_paths(simulations_dict: dict) -> None:
    """Raise FileNotFoundError early if any trajectory or topology is missing."""
    missing = []
    for sim_name, repeats in simulations_dict.items():
        for idx, (traj, top) in enumerate(repeats, start=1):
            for path, label in [(traj, "trajectory"), (top, "topology")]:
                if not os.path.exists(path):
                    missing.append(f"  [{sim_name} repeat {idx}] {label}: {path}")
    if missing:
        raise FileNotFoundError(
            "The following files were not found:\n" + "\n".join(missing)
        )


def _compute_dihedral_for_residue(args: tuple) -> str:
    """
    Compute dihedral angles for a single residue across all trajectory frames
    and save them to a CSV file.

    Parameters
    ----------
    args : tuple
        (sim_name, repeat_idx, top_path, traj_path, segid,
         resname, resid, atom_selector, dih_label)

    Returns
    -------
    str
        Status message prefixed with [DONE], [SKIP], or [FAIL].
    """
    sim_name, repeat_idx, top_path, traj_path, segid, resname, resid, atom_selector, dih_label = args

    if atom_selector is None:
        return f"[SKIP] {resname}{resid} – no {dih_label} definition"

    try:
        u = mda.Universe(top_path, traj_path)
        residue_atoms = u.select_atoms(f"segid {segid} and resid {resid}")
        dih_atoms = residue_atoms.select_atoms(atom_selector)

        if len(dih_atoms) < 4:
            return (
                f"[SKIP] {resname}{resid} – fewer than 4 atoms selected "
                f"for {dih_label} ({len(dih_atoms)} found)"
            )

        angles = []
        for _ts in u.trajectory:
            try:
                angle = np.round(
                    np.degrees(
                        mda.lib.distances.calc_dihedrals(
                            dih_atoms.positions[0],
                            dih_atoms.positions[1],
                            dih_atoms.positions[2],
                            dih_atoms.positions[3],
                        )
                    ),
                    3,
                )
                angles.append(angle)
            except Exception:
                continue  # skip problematic frames silently

        if not angles:
            return f"[SKIP] {resname}{resid} – no valid frames"

        filename = f"{sim_name}_repeat_{repeat_idx}_{dih_label}_resid_{resname}{resid}.csv"
        np.savetxt(filename, angles, fmt="%.3f", delimiter=",", comments="")
        return f"[DONE] {sim_name} repeat {repeat_idx} {dih_label} {resname}{resid}"

    except Exception as exc:
        return f"[FAIL] {sim_name} repeat {repeat_idx} {resname}{resid}: {exc}"


def _process_repeat(args: tuple) -> list:
    """
    Build the task list for one simulation repeat and dispatch to a worker pool.

    Parameters
    ----------
    args : tuple
        (sim_name, repeat_idx, traj_path, top_path, segid,
         dihedral_label, dihedral_definitions)

    Returns
    -------
    list of str
        Status messages for each residue task.
    """
    sim_name, repeat_idx, traj_path, top_path, segid, dihedral_label, dihedral_definitions = args

    try:
        u = mda.Universe(top_path, traj_path)
        residue_tasks = []

        for resname, atom_selector in dihedral_definitions.Resids.items():
            atoms = u.select_atoms(f"resname {resname} and segid {segid}")
            resids = {res.resid for res in atoms.residues}
            for resid in resids:
                residue_tasks.append((
                    sim_name, repeat_idx, top_path, traj_path,
                    segid, resname, resid, atom_selector, dihedral_label,
                ))

        if not residue_tasks:
            return [f"[SKIP] {sim_name} repeat {repeat_idx} – no matching residues in segid '{segid}'"]

        n_workers = max(1, cpu_count() // 2)
        with Pool(processes=n_workers) as pool:
            results = list(pool.imap_unordered(_compute_dihedral_for_residue, residue_tasks))
        return results

    except Exception as exc:
        return [f"[FAIL] {sim_name} repeat {repeat_idx}: {exc}"]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_parallel_extraction(
    simulations_dict: dict,
    segid: str,
    dihedral_label: str,
    dihedral_definitions,
) -> None:
    """
    Extract dihedral angles for all residues in all simulation repeats.

    Repeats are processed sequentially (with a progress bar); residues within
    each repeat are processed in parallel using half the available CPU cores.

    Parameters
    ----------
    simulations_dict : dict
        Mapping of simulation name → list of (trajectory_path, topology_path).
    segid : str
        Chain / segment identifier in the topology (e.g. "P0").
    dihedral_label : str
        "chi1" or "chi2".
    dihedral_definitions : class
        DihedralDefinitions_chi1 or DihedralDefinitions_chi2 from
        dihedral_definitions.py.
    """
    print(f"[INFO] Validating input file paths ...")
    _validate_paths(simulations_dict)

    repeat_args = []
    for sim_name, repeats in simulations_dict.items():
        for idx, (traj, top) in enumerate(repeats, start=1):
            repeat_args.append((sim_name, idx, traj, top, segid, dihedral_label, dihedral_definitions))

    n_repeats = len(repeat_args)
    print(f"[INFO] Processing {n_repeats} repeat(s) across {len(simulations_dict)} simulation(s).")

    all_messages = []
    for args in tqdm(repeat_args, desc="Extracting dihedrals"):
        messages = _process_repeat(args)
        all_messages.extend(messages)

    # Print summary counts
    done  = sum(1 for m in all_messages if m.startswith("[DONE]"))
    skip  = sum(1 for m in all_messages if m.startswith("[SKIP]"))
    fail  = sum(1 for m in all_messages if m.startswith("[FAIL]"))
    print(f"\n[EXTRACTION COMPLETE] DONE: {done}  SKIP: {skip}  FAIL: {fail}")

    if fail:
        print("\nFailed tasks:")
        for m in all_messages:
            if m.startswith("[FAIL]"):
                print(" ", m)
