---
title: AmberFlow - MD Simulation Pipeline
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
pinned: false
license: mit
app_port: 7860
---

# AmberFlow - Molecular Dynamics Simulation Pipeline

🧬 **AmberFlow** is a comprehensive web-based pipeline for preparing and setting up molecular dynamics (MD) simulations using the AMBER force field. This tool provides an intuitive interface for protein structure preparation, parameter generation, and simulation file creation.

## Features

### 🔬 Structure Preparation
- **Protein Loading**: Upload PDB files or fetch from RCSB PDB database
- **Structure Cleaning**: Remove water molecules, ions, and hydrogen atoms
- **Capping Groups**: Add ACE (N-terminal) and NME (C-terminal) capping groups
- **Ligand Handling**: Preserve and process ligands with automatic force field parameter generation
- **3D Visualization**: Interactive molecular viewer using NGL

### ⚙️ Simulation Parameters
- **Force Fields**: Support for ff14SB and ff19SB protein force fields
- **Water Models**: TIP3P and SPCE water models
- **System Setup**: Configurable box size and ion addition
- **Thermodynamics**: Temperature and pressure control

### 📋 Simulation Steps
- **Restrained Minimization**: Position-restrained energy minimization
- **Minimization**: Full system energy minimization
- **NPT Heating**: Temperature equilibration
- **NPT Equilibration**: Pressure and temperature equilibration
- **Production Run**: Configurable production MD simulation

### 📁 File Generation
- **AMBER Input Files**: Complete set of .in files for all simulation steps
- **Force Field Parameters**: Generated .prmtop and .inpcrd files
- **PBS Scripts**: HPC submission scripts
- **Analysis Scripts**: Post-simulation analysis tools

## Usage

1. **Load Protein Structure**
   - Upload a PDB file or enter a PDB ID to fetch from RCSB
   - View 3D structure and basic information

2. **Prepare Structure**
   - Configure structure preparation options
   - Remove unwanted components (water, ions, hydrogens)
   - Add capping groups for termini
   - Handle ligands if present

3. **Set Simulation Parameters**
   - Choose force field and water model
   - Configure system parameters
   - Set temperature and pressure

4. **Configure Simulation Steps**
   - Enable/disable simulation steps
   - Set step-specific parameters
   - Configure production run duration

5. **Generate Files**
   - Generate all simulation input files
   - Download files as ZIP archive
   - Preview generated files

## Technical Details

### Dependencies
- **MDAnalysis**: Structure manipulation and analysis
- **BioPython**: PDB file parsing
- **Flask**: Web framework
- **NGL Viewer**: 3D molecular visualization
- **AMBER Tools**: Force field parameter generation

### File Structure
```
AmberFlow/
├── app.py                 # Hugging Face Spaces entry point
├── requirements.txt       # Python dependencies
├── python/
│   ├── app.py            # Main Flask application
│   ├── structure_preparation.py
│   └── requirements.txt
├── html/
│   └── index.html        # Web interface
├── css/
│   └── styles.css        # Styling
├── js/
│   └── script.js         # Frontend logic
├── templates/            # AMBER input file templates
└── add_caps.py          # Capping group addition script
```

## Citation

If you use AmberFlow in your research, please cite:

```bibtex
@software{Amberflow2025,
  title={AmberFlow: Molecular Dynamics Simulation Pipeline},
  author={Hemant Nagar},
  year={2025},
  url={https://huggingface.co/spaces/hemantn/AmberFlow}
}
```

## Acknowledgments

- **Mohd Ibrahim** (Technical University of Munich) for the protein capping functionality (`add_caps.py`)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

- **Author**: Hemant Nagar
- **Email**: hn533621@ohio.edu

