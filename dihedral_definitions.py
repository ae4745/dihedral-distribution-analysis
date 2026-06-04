"""
dihedral_definitions.py
========================
Atom selections for χ1 and χ2 side-chain dihedral angles, covering the 20
standard amino acids plus common CHARMM/AMBER protonation-state variants.

Format
------
Each class exposes a ``Resids`` dict mapping a three-letter residue name (as
it appears in your topology) to an MDAnalysis atom-name selector string.
A value of ``None`` means the residue has no meaningful dihedral of that type
and will be silently skipped during extraction.

Adding a new residue type
-------------------------
1. Add an entry to ``DihedralDefinitions_chi1.Resids`` (and optionally chi2).
2. Use the four atom names that define the dihedral, space-separated, matching
   exactly the names in your topology file.
   Example:  'name N CA CB CG'
3. If the dihedral does not exist for that residue, set the value to ``None``.
"""


class DihedralDefinitions_chi1:
    """
    χ1 dihedral: N – CA – CB – (first side-chain heavy atom)

    Measures rotation about the CA–CB bond. Defined for all residues with a
    CB atom; GLY and ALA are included as backbone pseudo-dihedrals for
    completeness (N–CA–C–O) but are rarely informative.
    """

    Resids = {
        # Standard amino acids
        "ARG":  "name N CA CB CG",
        "ALA":  "name N CA C O",      # no true chi1; backbone pseudo-dihedral
        "ASN":  "name N CA CB CG",
        "ASP":  "name N CA CB CG",
        "CYS":  "name N CA CB SG",
        "GLN":  "name N CA CB CG",
        "GLU":  "name N CA CB CG",
        "GLY":  "name N CA C O",      # no true chi1; backbone pseudo-dihedral
        "HIS":  "name N CA CB CG",
        "ILE":  "name N CA CB CG1",
        "LEU":  "name N CA CB CG",
        "LYS":  "name N CA CB CG",
        "MET":  "name N CA CB CG",
        "PHE":  "name N CA CB CG",
        "PRO":  "name N CA CB CG",
        "SER":  "name N CA CB OG",
        "THR":  "name N CA CB OG1",
        "TRP":  "name N CA CB CG",
        "TYR":  "name N CA CB CG",
        "VAL":  "name N CA CB CG1",
        # CHARMM protonation-state variants
        "ASPP": "name N CA CB CG",    # protonated ASP
        "CYSP": "name N CA CB SG",    # protonated CYS
        "CYSG": "name N CA CB SG",    # CYS in disulfide
        "HSD":  "name N CA CB CG",    # CHARMM neutral HIS (δ)
        "HSE":  "name N CA CB CG",    # CHARMM neutral HIS (ε)
        "HSP":  "name N CA CB CG",    # CHARMM protonated HIS
        "SERD": "name N CA CB OG",    # deprotonated SER
        # AMBER protonation-state variants
        "HID":  "name N CA CB CG",    # AMBER neutral HIS (δ)
        "HIP":  "name N CA CB CG",    # AMBER protonated HIS
    }


class DihedralDefinitions_chi2:
    """
    χ2 dihedral: CA – CB – CG – (next heavy atom)

    Measures rotation about the CB–CG bond. Only defined for residues with at
    least three side-chain heavy atoms after CA. Residues with ``None`` are
    automatically skipped.
    """

    Resids = {
        # Standard amino acids
        "ARG":  "name CA CB CG CD",
        "ALA":  None,
        "ASN":  "name CA CB CG OD1",
        "ASP":  "name CA CB CG OD1",
        "CYS":  None,
        "GLN":  "name CA CB CG CD",
        "GLU":  "name CA CB CG CD",
        "GLY":  None,
        "HIS":  "name CA CB CG ND1",
        "ILE":  "name CA CB CG1 CD1",
        "LEU":  "name CA CB CG CD1",
        "LYS":  "name CA CB CG CD",
        "MET":  "name CA CB CG SD",
        "PHE":  "name CA CB CG CD1",
        "PRO":  "name CA CB CG CD",
        "SER":  None,
        "THR":  None,
        "TRP":  "name CA CB CG CD1",
        "TYR":  "name CA CB CG CD1",
        "VAL":  None,
        # CHARMM protonation-state variants
        "ASPP": "name CA CB CG OD1",
        "CYSP": None,
        "CYSG": None,
        "HSD":  "name CA CB CG ND1",
        "HSE":  "name CA CB CG ND1",
        "HSP":  "name CA CB CG ND1",
        "SERD": None,
        # AMBER protonation-state variants
        "HID":  "name CA CB CG ND1",
        "HIP":  "name CA CB CG ND1",
    }
