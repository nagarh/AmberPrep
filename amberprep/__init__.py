"""
AmberPrep: Web-based MD simulation pipeline with AMBER, ESMFold, docking, and PLUMED.

AmberPrep provides a complete workflow for:
- Protein structure loading and visualization
- Missing residue completion with ESMFold
- Structure preparation (cleaning, capping, chain/ligand selection)
- Ligand docking with AutoDock Vina + Meeko
- AMBER force field parameterization
- MD simulation file generation
- PLUMED collective variable configuration

Usage:
    # Run the web interface
    $ amberprep
    # or
    $ python -m amberprep

    # Import in Python
    from amberprep.app import app
    from amberprep.structure_preparation import prepare_structure

Requirements:
    - Python >= 3.10
    - Conda packages: ambertools, pymol-open-source, vina, openbabel, rdkit, gemmi
    - See README.md for full installation instructions

License: MIT
"""

__version__ = "0.1.3"
__author__ = "Hemant Nagar"
__email__ = "hn533621@ohio.edu"

# Expose key components for programmatic use
from amberprep.app import app

__all__ = ["app", "__version__"]
