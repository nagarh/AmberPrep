#!/usr/bin/env python3
"""
MD Simulation Pipeline - Flask Backend
Provides API endpoints for protein processing and file generation
"""

from flask import Flask, request, jsonify, send_file, render_template, send_from_directory, Response, stream_with_context
from flask_cors import CORS
import os
import json
import tempfile
import zipfile
from pathlib import Path
import requests
import subprocess
from Bio.PDB import PDBParser, PDBList
import logging
import html
from structure_preparation import prepare_structure, parse_structure_info
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from Fill_missing_residues import (
    get_pdb_id_from_pdb_file,
    detect_missing_residues,
    get_chain_sequences,
    run_esmfold,
    rebuild_pdb_with_esmfold,
    write_fasta_for_missing_chains,
    trim_residues_from_edges,
    trim_chains_sequences
)

app = Flask(__name__, 
            template_folder='../html',
            static_folder='../',
            static_url_path='')
CORS(app)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create output directory
OUTPUT_DIR = Path(__file__).parent.parent / "output"

def clean_and_create_output_folder():
    """Clean existing output folder and create a new one"""
    try:
        print(f"DEBUG: Starting cleanup. OUTPUT_DIR = {OUTPUT_DIR}")
        print(f"DEBUG: OUTPUT_DIR.exists() = {OUTPUT_DIR.exists()}")
        
        # Remove existing output folder if it exists
        if OUTPUT_DIR.exists():
            import shutil
            print(f"DEBUG: Removing existing output folder: {OUTPUT_DIR}")
            shutil.rmtree(OUTPUT_DIR)
            print(f"DEBUG: Successfully removed output folder")
            logger.info(f"Removed existing output folder: {OUTPUT_DIR}")
        
        # Create new output folder
        print(f"DEBUG: Creating new output folder: {OUTPUT_DIR}")
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        print(f"DEBUG: Successfully created output folder")
        logger.info(f"Created new output folder: {OUTPUT_DIR}")
        
        return True
    except Exception as e:
        print(f"DEBUG: Error in cleanup: {str(e)}")
        logger.error(f"Error cleaning output folder: {str(e)}")
        return False

class MDSimulationGenerator:
    """Handles MD simulation file generation and protein processing"""
    
    def __init__(self):
        self.pdb_parser = PDBParser(QUIET=True)
        self.pdb_list = PDBList()
    
    def fetch_pdb_structure(self, pdb_id):
        """Fetch PDB structure from RCSB"""
        try:
            # Download PDB file
            pdb_file = self.pdb_list.retrieve_pdb_file(pdb_id, pdir=OUTPUT_DIR, file_format='pdb')
            return str(pdb_file)
        except Exception as e:
            logger.error(f"Error fetching PDB {pdb_id}: {str(e)}")
            raise
    
    def parse_pdb_structure(self, pdb_file):
        """Parse PDB file and extract structure information"""
        try:
            structure = self.pdb_parser.get_structure('protein', pdb_file)
            
            # Extract basic information
            atom_count = 0
            chains = set()
            residues = set()
            
            for model in structure:
                for chain in model:
                    chains.add(chain.id)
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residues
                            residues.add(f"{residue.resname}{residue.id[1]}")
                        for atom in residue:
                            atom_count += 1
            
            return {
                'atom_count': atom_count,
                'chains': list(chains),
                'residue_count': len(residues),
                'structure_id': Path(pdb_file).stem.upper()
            }
        except Exception as e:
            logger.error(f"Error parsing PDB file: {str(e)}")
            raise
    
    def generate_mdp_file(self, params, step_type='production'):
        """Generate GROMACS MDP file for different simulation steps"""
        
        if step_type == 'restrained_min':
            return f"""; Restrained Minimization Parameters
integrator = steep
nsteps = {params['steps']['restrainedMin']['steps']}
emstep = 0.01
emtol = 1000

; Position restraints
define = -DPOSRES
refcoord_scaling = com

; Output control
nstxout = 100
nstenergy = 100
nstlog = 100

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = {params['cutoff']}

; Electrostatics
coulombtype = PME
rcoulomb = {params['cutoff']}
pme_order = {params['pmeOrder']}

; Van der Waals
vdwtype = Cut-off
rvdw = {params['cutoff']}
"""
        
        elif step_type == 'minimization':
            return f"""; Minimization Parameters
integrator = {params['steps']['minimization']['algorithm']}
nsteps = {params['steps']['minimization']['steps']}
emstep = 0.01
emtol = 1000

; Output control
nstxout = 100
nstenergy = 100
nstlog = 100

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = {params['cutoff']}

; Electrostatics
coulombtype = PME
rcoulomb = {params['cutoff']}
pme_order = {params['pmeOrder']}

; Van der Waals
vdwtype = Cut-off
rvdw = {params['cutoff']}
"""
        
        elif step_type == 'nvt':
            return f"""; NVT Equilibration Parameters
integrator = md
dt = {params['timestep']}
nsteps = {params['steps']['nvt']['steps']}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = {params['cutoff']}

; Electrostatics
coulombtype = PME
rcoulomb = {params['cutoff']}
pme_order = {params['pmeOrder']}

; Van der Waals
vdwtype = Cut-off
rvdw = {params['cutoff']}

; Temperature coupling
tcoupl = {params['couplingType']}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = {params['steps']['nvt']['temperature']} {params['steps']['nvt']['temperature']}

; Pressure coupling (disabled for NVT)
pcoupl = no

; Velocity generation
gen_vel = yes
gen_temp = {params['steps']['nvt']['temperature']}
gen_seed = -1
"""
        
        elif step_type == 'npt':
            return f"""; NPT Equilibration Parameters
integrator = md
dt = {params['timestep']}
nsteps = {params['steps']['npt']['steps']}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = {params['cutoff']}

; Electrostatics
coulombtype = PME
rcoulomb = {params['cutoff']}
pme_order = {params['pmeOrder']}

; Van der Waals
vdwtype = Cut-off
rvdw = {params['cutoff']}

; Temperature coupling
tcoupl = {params['couplingType']}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = {params['steps']['npt']['temperature']} {params['steps']['npt']['temperature']}

; Pressure coupling
pcoupl = {params['couplingType']}
pcoupltype = isotropic
tau_p = 2.0
ref_p = {params['steps']['npt']['pressure']}
compressibility = 4.5e-5

; Velocity generation
gen_vel = no
"""
        
        else:  # production
            return f"""; MD Simulation Parameters
; Generated by MD Simulation Pipeline

; Run parameters
integrator = md
dt = {params['timestep']}
nsteps = {params['steps']['production']['steps']}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = {params['cutoff']}

; Electrostatics
coulombtype = PME
rcoulomb = {params['cutoff']}
pme_order = {params['pmeOrder']}
fourierspacing = 0.16

; Van der Waals
vdwtype = Cut-off
rvdw = {params['cutoff']}

; Temperature coupling
tcoupl = {params['couplingType']}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = {params['temperature']} {params['temperature']}

; Pressure coupling
pcoupl = {params['couplingType']}
pcoupltype = isotropic
tau_p = 2.0
ref_p = {params['pressure']}
compressibility = 4.5e-5

; Dispersion correction
DispCorr = EnerPres

; Velocity generation
gen_vel = yes
gen_temp = {params['temperature']}
gen_seed = -1
"""
    
    def generate_pbs_script(self, protein_name, params):
        """Generate PBS script for HPC submission"""
        total_steps = params['steps']['production']['steps']
        time_in_ns = (total_steps * params['timestep']) / 1000
        
        return f"""#!/bin/bash
#PBS -N {protein_name}_md
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -q normal
#PBS -j oe

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Load required modules
module load gromacs/2023.2
module load intel/2021.4.0

# Set up environment
export OMP_NUM_THREADS=16
export GMX_MAXBACKUP=-1

# Simulation parameters
PROTEIN={protein_name}
STEPS={total_steps}
TIME_NS={time_in_ns:.2f}

echo "Starting MD simulation for $PROTEIN"
echo "Total simulation time: $TIME_NS ns"
echo "Job started at: $(date)"

# Run the simulation
./run_simulation.sh $PROTEIN

echo "Simulation completed at: $(date)"
echo "Results saved in output directory"
"""
    
    def generate_setup_script(self, protein_name, params):
        """Generate setup script for MD simulation"""
        return f"""#!/bin/bash
# Setup script for {protein_name} MD simulation
# Generated by MD Simulation Pipeline

set -e

PROTEIN={protein_name}
FORCE_FIELD={params['forceField']}
WATER_MODEL={params['waterModel']}

echo "Setting up MD simulation for $PROTEIN"

# Create output directory
mkdir -p output

# 1. Prepare protein structure
echo "Preparing protein structure..."
gmx pdb2gmx -f $PROTEIN.pdb -o $PROTEIN_processed.gro -p $PROTEIN.top -ff $FORCE_FIELD -water $WATER_MODEL

# 2. Define simulation box
echo "Defining simulation box..."
gmx editconf -f $PROTEIN_processed.gro -o $PROTEIN_box.gro -c -d {params['boxMargin']} -bt {params['boxType']}

# 3. Add solvent
echo "Adding solvent..."
gmx solvate -cp $PROTEIN_box.gro -cs spc216.gro -o $PROTEIN_solv.gro -p $PROTEIN.top

# 4. Add ions
echo "Adding ions..."
gmx grompp -f $PROTEIN_restrained.mdp -c $PROTEIN_solv.gro -p $PROTEIN.top -o $PROTEIN_ions.tpr
echo "SOL" | gmx genion -s $PROTEIN_ions.tpr -o $PROTEIN_final.gro -p $PROTEIN.top -pname NA -nname CL -neutral

echo "Setup completed successfully!"
echo "Ready to run simulation with: ./run_simulation.sh $PROTEIN"
"""
    
    def generate_analysis_script(self, protein_name):
        """Generate analysis script for MD simulation results"""
        return f"""#!/bin/bash
# Analysis script for {protein_name} MD simulation
# Generated by MD Simulation Pipeline

PROTEIN={protein_name}

echo "Analyzing MD simulation results for $PROTEIN"

# Create analysis directory
mkdir -p analysis

# 1. RMSD analysis
echo "Calculating RMSD..."
echo "Protein" | gmx rms -s $PROTEIN_final.tpr -f $PROTEIN_prod.xtc -o analysis/$PROTEIN_rmsd.xvg -tu ns

# 2. RMSF analysis
echo "Calculating RMSF..."
echo "Protein" | gmx rmsf -s $PROTEIN_final.tpr -f $PROTEIN_prod.xtc -o analysis/$PROTEIN_rmsf.xvg -res

# 3. Radius of gyration
echo "Calculating radius of gyration..."
echo "Protein" | gmx gyrate -s $PROTEIN_final.tpr -f $PROTEIN_prod.xtc -o analysis/$PROTEIN_gyrate.xvg

# 4. Hydrogen bonds
echo "Analyzing hydrogen bonds..."
echo "Protein" | gmx hbond -s $PROTEIN_final.tpr -f $PROTEIN_prod.xtc -num analysis/$PROTEIN_hbonds.xvg

# 5. Energy analysis
echo "Analyzing energies..."
gmx energy -f $PROTEIN_prod.edr -o analysis/$PROTEIN_energy.xvg

# 6. Generate plots
echo "Generating analysis plots..."
python3 plot_analysis.py $PROTEIN

echo "Analysis completed! Results saved in analysis/ directory"
"""

# Initialize the MD simulation generator
md_generator = MDSimulationGenerator()

@app.route('/api/fetch-pdb', methods=['POST'])
def fetch_pdb():
    """Fetch PDB structure from RCSB"""
    try:
        print("DEBUG: fetch-pdb endpoint called")
        data = request.get_json()
        pdb_id = data.get('pdb_id', '').upper()
        print(f"DEBUG: pdb_id = {pdb_id}")
        
        if not pdb_id or len(pdb_id) != 4:
            return jsonify({'error': 'Invalid PDB ID'}), 400
        
        # Clean and create new output folder for fresh start
        print("DEBUG: Calling clean_and_create_output_folder()")
        if not clean_and_create_output_folder():
            return jsonify({'error': 'Failed to clean output folder'}), 500
        print("DEBUG: Output folder cleanup completed successfully")
        
        # Fetch PDB structure
        pdb_file = md_generator.fetch_pdb_structure(pdb_id)
        
        # Parse structure information
        structure_info = md_generator.parse_pdb_structure(pdb_file)
        
        return jsonify({
            'success': True,
            'structure_info': structure_info,
            'pdb_file': pdb_file
        })
    
    except Exception as e:
        logger.error(f"Error fetching PDB: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/parse-pdb', methods=['POST'])
def parse_pdb():
    """Parse uploaded PDB file"""
    try:
        print("DEBUG: parse-pdb endpoint called")
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        print(f"DEBUG: Processing uploaded file: {file.filename}")
        
        # Clean and create new output folder for fresh start
        print("DEBUG: Calling clean_and_create_output_folder()")
        if not clean_and_create_output_folder():
            return jsonify({'error': 'Failed to clean output folder'}), 500
        print("DEBUG: Output folder cleanup completed successfully")
        
        # Save uploaded file temporarily
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        file.save(temp_file.name)
        
        # Parse structure information
        structure_info = md_generator.parse_pdb_structure(temp_file.name)
        
        # Clean up temporary file
        os.unlink(temp_file.name)
        
        return jsonify({
            'success': True,
            'structure_info': structure_info
        })
    
    except Exception as e:
        logger.error(f"Error parsing PDB: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/generate-files', methods=['POST'])
def generate_files():
    """Generate MD simulation files"""
    try:
        data = request.get_json()
        protein_name = data.get('protein_name', 'protein')
        simulation_params = data.get('simulation_params', {})
        
        # Generate all files
        files = {}
        
        # MDP files
        files[f'{protein_name}.mdp'] = md_generator.generate_mdp_file(simulation_params, 'production')
        files[f'{protein_name}_restrained.mdp'] = md_generator.generate_mdp_file(simulation_params, 'restrained_min')
        files[f'{protein_name}_min.mdp'] = md_generator.generate_mdp_file(simulation_params, 'minimization')
        files[f'{protein_name}_nvt.mdp'] = md_generator.generate_mdp_file(simulation_params, 'nvt')
        files[f'{protein_name}_npt.mdp'] = md_generator.generate_mdp_file(simulation_params, 'npt')
        files[f'{protein_name}_prod.mdp'] = md_generator.generate_mdp_file(simulation_params, 'production')
        
        # Scripts
        files[f'{protein_name}_simulation.pbs'] = md_generator.generate_pbs_script(protein_name, simulation_params)
        files[f'setup_{protein_name}.sh'] = md_generator.generate_setup_script(protein_name, simulation_params)
        files[f'analyze_{protein_name}.sh'] = md_generator.generate_analysis_script(protein_name)
        
        return jsonify({
            'success': True,
            'files': files
        })
    
    except Exception as e:
        logger.error(f"Error generating files: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/download-zip', methods=['POST'])
def download_zip():
    """Download all generated files as a ZIP archive"""
    try:
        data = request.get_json()
        files = data.get('files', {})
        
        # Create temporary ZIP file
        temp_zip = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')
        
        with zipfile.ZipFile(temp_zip.name, 'w') as zip_file:
            for filename, content in files.items():
                zip_file.writestr(filename, content)
        
        return send_file(
            temp_zip.name,
            as_attachment=True,
            download_name='md_simulation_files.zip',
            mimetype='application/zip'
        )
    
    except Exception as e:
        logger.error(f"Error creating ZIP file: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/get-solvated-protein', methods=['GET'])
def get_solvated_protein():
    """Get the solvated protein PDB file content"""
    try:
        solvated_file = os.path.join(OUTPUT_DIR, 'protein_solvated.pdb')
        
        if not os.path.exists(solvated_file):
            return jsonify({'success': False, 'error': 'Solvated protein file not found. Please generate files first.'})
        
        with open(solvated_file, 'r') as f:
            content = f.read()
        
        return jsonify({'success': True, 'content': content})
    except Exception as e:
        logger.error(f"Error reading solvated protein file: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/get-viewer-pdb', methods=['GET'])
def get_viewer_pdb():
    """Return a single PDB for viewer: start from protein_solvated.pdb and mark ligand residues as HETATM.
    Ligand residues are detected from 4_ligands_corrected*.pdb files by (resname, chain, resi) tuples; if chains/resi not present, fallback to resname matching.
    """
    try:
        solvated_path = OUTPUT_DIR / 'protein_solvated.pdb'
        # Find all corrected ligand files (support multiple ligands)
        # Exclude OpenBabel output files (4_ligands_corrected_obabel_*.pdb)
        lig_paths = sorted([f for f in OUTPUT_DIR.glob('4_ligands_corrected_*.pdb') if "_obabel_" not in f.name])
        # Fallback to single file for backward compatibility
        if not lig_paths:
            single_lig_path = OUTPUT_DIR / '4_ligands_corrected.pdb'
            if single_lig_path.exists():
                lig_paths = [single_lig_path]
        viewer_out = OUTPUT_DIR / 'viewer_protein_with_ligand.pdb'

        if not solvated_path.exists():
            return jsonify({'success': False, 'error': 'protein_solvated.pdb not found'}), 400

        # Build ligand index from all corrected ligand PDB files if present
        ligand_keys = set()
        ligand_resnames = set()
        for lig_path in lig_paths:
            if lig_path.exists():
                with open(lig_path, 'r') as lf:
                    for line in lf:
                        if line.startswith(('ATOM', 'HETATM')):
                            resn = line[17:20].strip()
                            chain = line[21:22].strip()
                            resi = line[22:26].strip()
                            ligand_resnames.add(resn)
                            if chain and resi:
                                ligand_keys.add((resn, chain, resi))

        # Rewrite solvated file marking matching ligand residues and ions (NA/CL) as HETATM
        out_lines = []
        with open(solvated_path, 'r') as sf:
            for line in sf:
                if line.startswith(('ATOM', 'HETATM')):
                    resn = line[17:20].strip()
                    chain = line[21:22].strip()
                    resi = line[22:26].strip()
                    is_match = False
                    is_ion = resn in { 'NA', 'CL' }
                    if (resn, chain, resi) in ligand_keys:
                        is_match = True
                    elif resn in ligand_resnames:
                        # Fallback by residue name only
                        is_match = True
                    if is_match or is_ion:
                        # Force to HETATM
                        out_lines.append('HETATM' + line[6:])
                    else:
                        out_lines.append(line)
                else:
                    out_lines.append(line)

        # Save combined viewer file (optional but useful for debugging)
        try:
            with open(viewer_out, 'w') as vf:
                vf.writelines(out_lines)
        except Exception:
            pass

        return jsonify({'success': True, 'content': ''.join(out_lines)})
    except Exception as e:
        logger.error(f"Error generating viewer PDB: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/view-pdb')
def view_pdb_html():
    """Serve PDB file as HTML page for instant viewing"""
    try:
        viewer_out = OUTPUT_DIR / 'viewer_protein_with_ligand.pdb'
        solvated_path = OUTPUT_DIR / 'protein_solvated.pdb'
        # Find all corrected ligand files (support multiple ligands)
        # Exclude OpenBabel output files (4_ligands_corrected_obabel_*.pdb)
        lig_paths = sorted([f for f in OUTPUT_DIR.glob('4_ligands_corrected_*.pdb') if "_obabel_" not in f.name])
        # Fallback to single file for backward compatibility
        if not lig_paths:
            single_lig_path = OUTPUT_DIR / '4_ligands_corrected.pdb'
            if single_lig_path.exists():
                lig_paths = [single_lig_path]
        
        # If viewer file doesn't exist, generate it first
        if not viewer_out.exists():
            if not solvated_path.exists():
                return f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <title>Error - PDB Not Found</title>
                    <style>
                        body {{ font-family: Arial, sans-serif; padding: 40px; text-align: center; }}
                        .error {{ color: #dc3545; font-size: 18px; }}
                    </style>
                </head>
                <body>
                    <div class="error">
                        <h1>PDB File Not Found</h1>
                        <p>Please complete the structure preparation steps first.</p>
                    </div>
                </body>
                </html>
                """, 404
            
            # Generate the file directly (same logic as get_viewer_pdb but without JSON response)
            try:
                # Build ligand index from all corrected ligand PDB files if present
                ligand_keys = set()
                ligand_resnames = set()
                for lig_path in lig_paths:
                    if lig_path.exists():
                        with open(lig_path, 'r') as lf:
                            for line in lf:
                                if line.startswith(('ATOM', 'HETATM')):
                                    resn = line[17:20].strip()
                                    chain = line[21:22].strip()
                                    resi = line[22:26].strip()
                                    ligand_resnames.add(resn)
                                    if chain and resi:
                                        ligand_keys.add((resn, chain, resi))

                # Rewrite solvated file marking matching ligand residues and ions (NA/CL) as HETATM
                out_lines = []
                with open(solvated_path, 'r') as sf:
                    for line in sf:
                        if line.startswith(('ATOM', 'HETATM')):
                            resn = line[17:20].strip()
                            chain = line[21:22].strip()
                            resi = line[22:26].strip()
                            is_match = False
                            is_ion = resn in { 'NA', 'CL' }
                            if (resn, chain, resi) in ligand_keys:
                                is_match = True
                            elif resn in ligand_resnames:
                                # Fallback by residue name only
                                is_match = True
                            if is_match or is_ion:
                                # Force to HETATM
                                out_lines.append('HETATM' + line[6:])
                            else:
                                out_lines.append(line)
                        else:
                            out_lines.append(line)

                # Save combined viewer file
                with open(viewer_out, 'w') as vf:
                    vf.writelines(out_lines)
            except Exception as e:
                logger.error(f"Error generating viewer PDB: {str(e)}")
                return f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <title>Error</title>
                    <style>
                        body {{ font-family: Arial, sans-serif; padding: 40px; text-align: center; }}
                        .error {{ color: #dc3545; font-size: 18px; }}
                    </style>
                </head>
                <body>
                    <div class="error">
                        <h1>Error Generating PDB</h1>
                        <p>Could not generate viewer PDB file: {html.escape(str(e))}</p>
                    </div>
                </body>
                </html>
                """, 500
        
        # Read PDB content
        with open(viewer_out, 'r') as f:
            pdb_content = f.read()
        
        # Escape HTML special characters
        escaped_content = html.escape(pdb_content)
        
        # Create HTML page
        html_page = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Viewer PDB File</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.4;
            background: #f8f9fa;
            padding: 20px;
        }}
        .header {{
            background: white;
            padding: 15px 20px;
            margin-bottom: 15px;
            border-radius: 4px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .header h1 {{
            font-size: 18px;
            color: #333;
        }}
        .pdb-content {{
            background: white;
            padding: 20px;
            border-radius: 4px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow-x: auto;
            white-space: pre;
            word-wrap: normal;
        }}
        .info {{
            color: #666;
            font-size: 11px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>📄 Viewer PDB File</h1>
        <div class="info">File: viewer_protein_with_ligand.pdb</div>
    </div>
    <div class="pdb-content">{escaped_content}</div>
</body>
</html>"""
        
        return html_page, 200, {'Content-Type': 'text/html; charset=utf-8'}
    except Exception as e:
        logger.error(f"Error serving PDB as HTML: {str(e)}")
        return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Error</title>
            <style>
                body {{ font-family: Arial, sans-serif; padding: 40px; text-align: center; }}
                .error {{ color: #dc3545; font-size: 18px; }}
            </style>
        </head>
        <body>
            <div class="error">
                <h1>Error Loading PDB</h1>
                <p>{html.escape(str(e))}</p>
            </div>
        </body>
        </html>
        """, 500

@app.route('/api/get-corrected-ligands', methods=['GET'])
def get_corrected_ligands():
    """Get the corrected ligand PDB file content if present (combines all ligands)"""
    try:
        # Find all corrected ligand files (support multiple ligands)
        # Exclude OpenBabel output files (4_ligands_corrected_obabel_*.pdb)
        ligand_files = sorted([f for f in OUTPUT_DIR.glob('4_ligands_corrected_*.pdb') if "_obabel_" not in f.name])
        # Fallback to single file for backward compatibility
        if not ligand_files:
            single_lig_file = OUTPUT_DIR / '4_ligands_corrected.pdb'
            if single_lig_file.exists():
                ligand_files = [single_lig_file]
        
        if not ligand_files:
            # Return success with exists flag false so frontend can decide gracefully
            return jsonify({'success': True, 'exists': False, 'content': ''})
        
        # Read and normalize records to HETATM for viewer compatibility, combine all ligands
        normalized_lines = []
        for ligand_file in ligand_files:
            with open(ligand_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        # Replace record name to HETATM, preserve fixed-width columns
                        normalized_lines.append('HETATM' + line[6:])
                    elif line.startswith('HETATM'):
                        normalized_lines.append(line)
                    elif line.strip() == 'END' and ligand_file != ligand_files[-1]:
                        # Skip END for intermediate ligands, keep only for last
                        continue
                    elif line.strip() and not line.startswith(('CRYST', 'REMARK', 'HEADER')):
                        normalized_lines.append(line)
        
        # Ensure we have an END at the end
        if normalized_lines and not normalized_lines[-1].strip() == 'END':
            normalized_lines.append('END\n')
        
        content = ''.join(normalized_lines)
        return jsonify({'success': True, 'exists': True, 'content': content})
    except Exception as e:
        logger.error(f"Error reading corrected ligand file: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/get-aligned-ligands', methods=['GET'])
def get_aligned_ligands():
    """Return ligand coordinates aligned to protein_solvated.pdb frame using PyMOL transforms."""
    try:
        solvated_file = OUTPUT_DIR / 'protein_solvated.pdb'
        tleap_ready = OUTPUT_DIR / 'tleap_ready.pdb'
        # Find all corrected ligand files (support multiple ligands)
        # Exclude OpenBabel output files (4_ligands_corrected_obabel_*.pdb)
        ligand_files = sorted([f for f in OUTPUT_DIR.glob('4_ligands_corrected_*.pdb') if "_obabel_" not in f.name])
        # Fallback to single file for backward compatibility
        if not ligand_files:
            single_lig_file = OUTPUT_DIR / '4_ligands_corrected.pdb'
            if single_lig_file.exists():
                ligand_files = [single_lig_file]

        if not solvated_file.exists():
            return jsonify({'success': False, 'error': 'protein_solvated.pdb not found'}), 400
        if not tleap_ready.exists():
            return jsonify({'success': False, 'error': 'tleap_ready.pdb not found'}), 400
        if not ligand_files:
            return jsonify({'success': True, 'exists': False, 'content': ''})
        
        # Use first ligand file for PyMOL alignment (or combine them if needed)
        ligand_file = ligand_files[0]

        # Create temp output path
        aligned_lig = OUTPUT_DIR / 'ligand_aligned_for_preview.pdb'
        try:
            if aligned_lig.exists():
                aligned_lig.unlink()
        except Exception:
            pass

        # PyMOL script: load solvated, load tlready (protein+lig), align tlready protein to solvated protein, then save transformed ligand
        pymol_script = f"""
import pymol
pymol.finish_launching(['pymol','-qc'])
from pymol import cmd
cmd.load('{solvated_file.as_posix()}', 'solv')
cmd.load('{tleap_ready.as_posix()}', 'prep')
cmd.load('{ligand_file.as_posix()}', 'lig')
# Align prepared protein to solvated protein; use CA atoms to be robust
cmd.align('prep and polymer.protein and name CA', 'solv and polymer.protein and name CA')
# Apply same transform implicitly affects 'prep' object; we saved ligand as separate object, so match matrices
mat = cmd.get_object_matrix('prep')
cmd.set_object_matrix('lig', mat)
# Save ligand in aligned frame, as HETATM
cmd.alter('lig', 'type="HETATM"')
cmd.save('{aligned_lig.as_posix()}', 'lig')
cmd.quit()
"""

        # Run PyMOL inline
        result = subprocess.run(['python3', '-c', pymol_script], capture_output=True, text=True, cwd=str(OUTPUT_DIR))
        if result.returncode != 0:
            return jsonify({'success': False, 'error': f'PyMOL alignment failed: {result.stderr}'}), 500

        if not aligned_lig.exists():
            return jsonify({'success': False, 'error': 'Aligned ligand file was not produced'}), 500

        # Read and return content
        normalized_lines = []
        with open(aligned_lig, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    normalized_lines.append('HETATM' + line[6:])
                else:
                    normalized_lines.append(line)
        content = ''.join(normalized_lines)
        return jsonify({'success': True, 'exists': True, 'content': content})
    except Exception as e:
        logger.error(f"Error aligning ligands: {str(e)}")
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/viewer/<filename>')
def viewer(filename):
    """Serve NGL viewer page"""
    # Check if the file exists, if not, try to generate it
    file_path = OUTPUT_DIR / filename
    if not file_path.exists():
        # Try to generate the viewer PDB if it's the specific file we need
        if filename == 'viewer_protein_with_ligand.pdb':
            try:
                # Call the get_viewer_pdb function to generate the file
                result = get_viewer_pdb()
                if result[1] == 200:  # Success
                    pass  # File should now exist
            except:
                pass  # Continue anyway
    
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>NGL Viewer - {filename}</title>
        <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js"></script>
        <style>
            body {{ margin: 0; padding: 0; font-family: Arial, sans-serif; }}
            #viewport {{ width: 100%; height: 100vh; }}
            .header {{ background: #f8f9fa; padding: 10px; border-bottom: 1px solid #ddd; }}
            .controls {{ padding: 10px; background: #f8f9fa; }}
            .btn {{ padding: 8px 16px; margin: 5px; border: none; border-radius: 4px; cursor: pointer; }}
            .btn-primary {{ background: #007bff; color: white; }}
            .btn-secondary {{ background: #6c757d; color: white; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h3>🧬 3D Structure Viewer - {filename}</h3>
        </div>
        <div id="viewport"></div>
        <div class="controls">
            <button class="btn btn-primary" onclick="resetView()">Reset View</button>
            <button class="btn btn-secondary" onclick="toggleRepresentation()">Toggle Style</button>
            <button class="btn btn-secondary" onclick="toggleSpin()">Toggle Spin</button>
        </div>
        <script>
            let stage;
            let currentRepresentation = 'cartoon';
            let isSpinning = false;

            async function initViewer() {{
                try {{
                    // Check if file exists first
                    const response = await fetch("/output/{filename}");
                    if (!response.ok) {{
                        throw new Error(`File not found: ${{response.status}} ${{response.statusText}}`);
                    }}
                    
                    stage = new NGL.Stage("viewport", {{ backgroundColor: "white" }});
                    
                    const component = await stage.loadFile("/output/{filename}");
                    
                    // Add cartoon representation for protein
                    component.addRepresentation("cartoon", {{
                        sele: "protein",
                        colorScheme: "chainname",
                        opacity: 0.9
                    }});

                    // Add ball and stick for water molecules
                    component.addRepresentation("ball+stick", {{
                        sele: "water",
                        color: "cyan",
                        colorScheme: "uniform",
                        radius: 0.1
                    }});

                    // Add ball and stick for ligands
                    component.addRepresentation("ball+stick", {{
                        sele: "hetero",
                        color: "element",
                        radius: 0.15
                    }});

                    stage.autoView();
                }} catch (error) {{
                    console.error('Error loading structure:', error);
                    document.getElementById('viewport').innerHTML = 
                        '<div style="padding: 50px; text-align: center; color: #dc3545;">' +
                        '<h3>Error loading structure</h3><p>' + error.message + '</p>' +
                        '<p>Make sure the file exists in the output directory.</p></div>';
                }}
            }}

            function resetView() {{
                if (stage) stage.autoView();
            }}

            function toggleRepresentation() {{
                if (!stage) return;
                const components = stage.compList;
                if (components.length === 0) return;

                const component = components[0];
                component.removeAllRepresentations();

                if (currentRepresentation === 'cartoon') {{
                    component.addRepresentation("ball+stick", {{
                        color: "element",
                        radius: 0.15
                    }});
                    currentRepresentation = 'ball+stick';
                }} else {{
                    component.addRepresentation("cartoon", {{
                        sele: "protein",
                        colorScheme: "chainname",
                        opacity: 0.9
                    }});
                    component.addRepresentation("ball+stick", {{
                        sele: "water",
                        color: "cyan",
                        colorScheme: "uniform",
                        radius: 0.1
                    }});
                    component.addRepresentation("ball+stick", {{
                        sele: "hetero",
                        color: "element",
                        radius: 0.15
                    }});
                    currentRepresentation = 'cartoon';
                }}
            }}

            function toggleSpin() {{
                if (!stage) return;
                isSpinning = !isSpinning;
                stage.setSpin(isSpinning);
            }}

            // Initialize when page loads
            document.addEventListener('DOMContentLoaded', initViewer);
        </script>
    </body>
    </html>
    """

@app.route('/output/<path:filename>')
def serve_output(filename):
    """Serve output files"""
    # Debug: print available files
    print(f"Requested file: {filename}")
    print(f"Full path: {OUTPUT_DIR / filename}")
    print(f"File exists: {(OUTPUT_DIR / filename).exists()}")
    print(f"Files in output dir: {list(OUTPUT_DIR.iterdir()) if OUTPUT_DIR.exists() else 'Directory not found'}")
    
    if not (OUTPUT_DIR / filename).exists():
        abort(404)
    
    return send_from_directory(OUTPUT_DIR, filename)

@app.route('/')
def index():
    """Serve the main HTML page"""
    return render_template('index.html')

@app.route('/<path:filename>')
def serve_static(filename):
    """Serve static files (CSS, JS, etc.)"""
    return send_from_directory('../', filename)

@app.route('/api/prepare-structure', methods=['POST'])
def prepare_structure_endpoint():
    """Prepare protein structure for AMBER"""
    try:
        data = request.get_json()
        pdb_content = data.get('pdb_content', '')
        options = data.get('options', {})
        
        # Check if completed structure exists and use it instead
        complete_structure_path = OUTPUT_DIR / "0_complete_structure.pdb"
        if complete_structure_path.exists():
            logger.info("Using completed structure (0_complete_structure.pdb) for preparation")
            with open(complete_structure_path, 'r') as f:
                pdb_content = f.read()
        elif not pdb_content:
            return jsonify({'error': 'No PDB content provided and no completed structure found'}), 400
        
        # Prepare structure
        result = prepare_structure(pdb_content, options)
        
        return jsonify({
            'success': True,
            'prepared_structure': result['prepared_structure'],
            'original_atoms': result['original_atoms'],
            'prepared_atoms': result['prepared_atoms'],
            'removed_components': result['removed_components'],
            'added_capping': result['added_capping'],
            'preserved_ligands': result['preserved_ligands'],
            'ligand_present': result.get('ligand_present', False),
            'separate_ligands': result.get('separate_ligands', False),
            'ligand_content': result.get('ligand_content', '')
        })
    
    except Exception as e:
        logger.error(f"Error preparing structure: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/parse-structure', methods=['POST'])
def parse_structure_endpoint():
    """Parse structure information"""
    try:
        data = request.get_json()
        pdb_content = data.get('pdb_content', '')
        
        if not pdb_content:
            return jsonify({'error': 'No PDB content provided'}), 400
        
        # Parse structure
        structure_info = parse_structure_info(pdb_content)
        
        return jsonify({
            'success': True,
            'structure_info': structure_info
        })
    
    except Exception as e:
        logger.error(f"Error parsing structure: {str(e)}")
        return jsonify({'error': str(e)}), 500

def _format_log(message, log_type='info'):
    """Helper function to format log message for SSE"""
    data = json.dumps({'type': log_type, 'message': message})
    return f"data: {data}\n\n"

@app.route('/api/generate-ligand-ff', methods=['POST'])
@stream_with_context
def generate_ligand_ff():
    """Generate force field parameters for multiple ligands with streaming logs"""
    def generate():
        try:
            data = request.get_json()
            force_field = data.get('force_field', 'gaff2')
            
            # Determine the s parameter based on force field
            s_param = 2 if force_field == 'gaff2' else 1
            
            yield _format_log(f"Working directory: {os.getcwd()}")
            yield _format_log(f"Output directory: {OUTPUT_DIR}")
            
            # Find all individual ligand files (4_ligands_corrected_1.pdb, 4_ligands_corrected_2.pdb, etc.)
            # Exclude OpenBabel output files (4_ligands_corrected_obabel_*.pdb)
            ligand_files = sorted([f for f in OUTPUT_DIR.glob("4_ligands_corrected_*.pdb") if "_obabel_" not in f.name])
            
            if not ligand_files:
                # Fallback: check for single ligand file (backward compatibility)
                single_ligand_pdb = OUTPUT_DIR / "4_ligands_corrected.pdb"
                if single_ligand_pdb.exists():
                    ligand_files = [single_ligand_pdb]
                else:
                    yield _format_log('Ligand PDB file(s) not found. Please prepare structure with ligands first.', 'error')
                    yield f"data: {json.dumps({'type': 'complete', 'success': False, 'error': 'Ligand PDB file(s) not found. Please prepare structure with ligands first.'})}\n\n"
                    return
            
            yield _format_log(f"Found {len(ligand_files)} ligand file(s) to process")
            
            import re
            processed_ligands = []
            errors = []
            
            # Step 1: Extract residue names and group ligands by residue name
            ligand_by_resname = {}  # Maps residue name to list of (ligand_pdb, ligand_num) tuples
            resname_to_ligand_num = {}  # Maps residue name to the ligand_num we'll use for processing
            
            for i, ligand_pdb in enumerate(ligand_files, 1):
                ligand_num = i
                # Extract number from filename if available (e.g., 4_ligands_corrected_1.pdb -> 1)
                match = re.search(r'_(\d+)\.pdb$', ligand_pdb.name)
                if match:
                    ligand_num = int(match.group(1))
                
                # Extract residue name from this ligand file
                resname = get_residue_name_from_pdb(ligand_pdb)
                if not resname:
                    yield _format_log(f"Warning: Could not extract residue name from {ligand_pdb.name}, using LIG{ligand_num}", 'warning')
                    resname = f"LIG{ligand_num}"
                
                # Group by residue name
                if resname not in ligand_by_resname:
                    ligand_by_resname[resname] = []
                    resname_to_ligand_num[resname] = ligand_num  # Use first occurrence's number
                ligand_by_resname[resname].append((ligand_pdb, ligand_num))
            
            yield _format_log(f"Found {len(ligand_by_resname)} unique ligand residue name(s): {', '.join(sorted(ligand_by_resname.keys()))}")
            
            # Step 2: Process each unique residue name only once
            for resname, ligand_list in ligand_by_resname.items():
                # Use the first ligand file for this residue name
                ligand_pdb, ligand_num = ligand_list[0]
                
                # If there are multiple occurrences, log it
                if len(ligand_list) > 1:
                    other_nums = [num for _, num in ligand_list[1:]]
                    yield _format_log(f"Residue {resname} appears {len(ligand_list)} times (ligand files: {ligand_num}, {', '.join(map(str, other_nums))})", 'info')
                    yield _format_log(f"Processing {resname} once using ligand file {ligand_num}, skipping duplicates", 'info')
                
                # Use residue name for output files to avoid conflicts
                ligand_mol2 = OUTPUT_DIR / f"{resname}.mol2"
                ligand_frcmod = OUTPUT_DIR / f"{resname}.frcmod"
                
                yield _format_log(f"\n{'='*60}")
                yield _format_log(f"Processing ligand {resname} (from file {ligand_pdb.name})")
                yield _format_log(f"{'='*60}")
                
                # Step 1: Calculate net charge using awk
                yield _format_log(f"Step 1: Calculating net charge for ligand {resname}...")
                awk_cmd = "awk '/^HETATM/ {if($NF ~ /[A-Z][0-9]-$/) charge--; if($NF ~ /[A-Z][0-9]\\+$/) charge++} END {print \"Net charge:\", charge+0}'"
                cmd1 = f"{awk_cmd} {ligand_pdb}"
                
                try:
                    result = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
                    output = result.stdout.strip()
                    yield _format_log(f"Awk output: '{output}'")
                    
                    net_charge_match = re.search(r'Net charge:\s*(-?\d+)', output)
                    if net_charge_match:
                        net_charge = int(net_charge_match.group(1))
                        yield _format_log(f"Calculated net charge: {net_charge}")
                    else:
                        yield _format_log("Could not extract net charge from awk output, using 0", 'warning')
                        net_charge = 0
                except Exception as e:
                    yield _format_log(f"Error running awk command: {e}, using net charge 0", 'error')
                    net_charge = 0
                
                # Step 2: Run antechamber with streaming output
                yield _format_log(f"Step 2: Running antechamber for ligand {resname} with net charge {net_charge}...")
                cmd2 = f"antechamber -i {ligand_pdb.name} -fi pdb -o {ligand_mol2.name} -fo mol2 -c bcc -at {force_field} -nc {net_charge}"
                yield _format_log(f"Running command: {cmd2}")
                
                # Stream antechamber output in real-time
                process = subprocess.Popen(cmd2, shell=True, cwd=str(OUTPUT_DIR), 
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                                         text=True, bufsize=1, universal_newlines=True)
                
                for line in iter(process.stdout.readline, ''):
                    if line:
                        yield _format_log(line.strip())
                
                process.wait()
                return_code = process.returncode
                
                yield _format_log(f"antechamber return code: {return_code}")
                
                if return_code != 0:
                    error_msg = f'antechamber failed for ligand {resname} with net charge {net_charge}'
                    yield _format_log(f"ERROR: {error_msg}", 'error')
                    errors.append(error_msg)
                    continue
                
                # Step 3: Run parmchk2 with streaming output
                yield _format_log(f"Step 3: Running parmchk2 for ligand {resname}...")
                cmd3 = f"parmchk2 -i {ligand_mol2.name} -f mol2 -o {ligand_frcmod.name} -a Y -s {s_param}"
                yield _format_log(f"Running command: {cmd3}")
                
                # Stream parmchk2 output in real-time
                process = subprocess.Popen(cmd3, shell=True, cwd=str(OUTPUT_DIR), 
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                                         text=True, bufsize=1, universal_newlines=True)
                
                for line in iter(process.stdout.readline, ''):
                    if line:
                        yield _format_log(line.strip())
                
                process.wait()
                return_code = process.returncode
                
                yield _format_log(f"parmchk2 return code: {return_code}")
                
                if return_code != 0:
                    error_msg = f'parmchk2 failed for ligand {resname}'
                    yield _format_log(f"ERROR: {error_msg}", 'error')
                    errors.append(error_msg)
                    continue
                
                # Check if files were generated successfully
                if ligand_mol2.exists() and ligand_frcmod.exists():
                    processed_ligands.append({
                        'resname': resname,
                        'ligand_num': ligand_num,
                        'net_charge': net_charge,
                        'files': {
                            'pdb': str(ligand_pdb),
                            'mol2': str(ligand_mol2),
                            'frcmod': str(ligand_frcmod)
                        },
                        'duplicate_files': [str(pdb) for pdb, num in ligand_list[1:]] if len(ligand_list) > 1 else []
                    })
                    yield _format_log(f"✅ Successfully processed ligand {resname}", 'success')
                else:
                    error_msg = f'Force field generation failed for ligand {resname} - output files not created'
                    yield _format_log(f"ERROR: {error_msg}", 'error')
                    errors.append(error_msg)
            
            if not processed_ligands:
                error_msg = f'Failed to process any ligands. Errors: {"; ".join(errors)}'
                yield _format_log(error_msg, 'error')
                yield f"data: {json.dumps({'type': 'complete', 'success': False, 'error': error_msg})}\n\n"
                return
            
            # Send final result
            result_data = {
                'type': 'complete',
                'success': True,
                'message': f'Successfully processed {len(processed_ligands)} ligand(s) with force field {force_field}',
                'ligands': processed_ligands,
                'errors': errors if errors else None
            }
            yield f"data: {json.dumps(result_data)}\n\n"
            
        except Exception as e:
            logger.error(f"Error generating ligand force field: {str(e)}")
            yield _format_log(f'Internal server error: {str(e)}', 'error')
            yield f"data: {json.dumps({'type': 'complete', 'success': False, 'error': f'Internal server error: {str(e)}'})}\n\n"
    
    return Response(generate(), mimetype='text/event-stream')

@app.route('/api/calculate-net-charge', methods=['POST'])
def calculate_net_charge():
    """Calculate net charge of the system using tleap"""
    try:
        # Check if structure is prepared
        tleap_ready_file = OUTPUT_DIR / "tleap_ready.pdb"
        if not tleap_ready_file.exists():
            return jsonify({'error': 'Structure not prepared. Please prepare structure first.'}), 400
        
        # Check if ligands are present - look for residue-named files first, then fallback to numbered files
        ligand_mol2_files = []
        ligand_frcmod_files = []
        ligand_resname_map = {}  # Maps residue name to (mol2_file, frcmod_file)
        
        # First, try to find residue-named files (e.g., O9C.mol2, O9C.frcmod)
        unique_resnames = get_all_ligand_residue_names()
        for resname in unique_resnames:
            mol2_file = OUTPUT_DIR / f"{resname}.mol2"
            frcmod_file = OUTPUT_DIR / f"{resname}.frcmod"
            if mol2_file.exists() and frcmod_file.exists():
                ligand_resname_map[resname] = (mol2_file, frcmod_file)
                ligand_mol2_files.append(mol2_file)
                ligand_frcmod_files.append(frcmod_file)
        
        # Fallback: check for numbered files (backward compatibility)
        if not ligand_mol2_files:
            numbered_mol2 = sorted(OUTPUT_DIR.glob("4_ligands_corrected_*.mol2"))
            numbered_frcmod = sorted(OUTPUT_DIR.glob("4_ligands_corrected_*.frcmod"))
            if numbered_mol2 and numbered_frcmod:
                ligand_mol2_files = numbered_mol2
                ligand_frcmod_files = numbered_frcmod
                # Try to map to residue names
                resnames = get_all_ligand_residue_names()
                for i, (mol2_file, frcmod_file) in enumerate(zip(ligand_mol2_files, ligand_frcmod_files)):
                    # Extract residue name from mol2 file if possible
                    resname = get_residue_name_from_mol2(mol2_file) if mol2_file.exists() else None
                    if not resname:
                        # Try to get from tleap_ready.pdb
                        if resnames and i < len(resnames):
                            resname = resnames[i]
                        else:
                            resname = f"LIG{len(ligand_resname_map) + 1}"
                    # Only add if not already in map (avoid duplicates)
                    if resname not in ligand_resname_map:
                        ligand_resname_map[resname] = (mol2_file, frcmod_file)
        
        # Final fallback: single ligand file (backward compatibility)
        if not ligand_mol2_files:
            single_mol2 = OUTPUT_DIR / "4_ligands_corrected.mol2"
            single_frcmod = OUTPUT_DIR / "4_ligands_corrected.frcmod"
            if single_mol2.exists() and single_frcmod.exists():
                ligand_mol2_files = [single_mol2]
                ligand_frcmod_files = [single_frcmod]
                resname = get_all_ligand_residue_names()
                if resname:
                    ligand_resname_map[resname[0]] = (single_mol2, single_frcmod)
                else:
                    ligand_resname_map["LIG"] = (single_mol2, single_frcmod)
        
        ligand_present = len(ligand_mol2_files) > 0 and len(ligand_frcmod_files) > 0
        
        # Create dynamic tleap input file
        tleap_input = OUTPUT_DIR / "calc_charge_on_system.in"
        
        # Get the selected force field from the request
        data = request.get_json() if request.get_json() else {}
        selected_force_field = data.get('force_field', 'ff14SB')
        
        with open(tleap_input, 'w') as f:
            f.write(f"source leaprc.protein.{selected_force_field}\n")
            f.write("source leaprc.gaff2\n\n")
            
            if ligand_present:
                # Load each unique ligand parameter and structure only once
                # Use sorted to ensure consistent ordering
                for resname in sorted(ligand_resname_map.keys()):
                    mol2_file, frcmod_file = ligand_resname_map[resname]
                    f.write(f"loadamberparams {frcmod_file.name}\n")
                    f.write(f"{resname} = loadmol2 {mol2_file.name}\n")
                f.write("\n")
            
            f.write("x = loadpdb tleap_ready.pdb\n\n")
            f.write("charge x\n\n")
            f.write("quit\n")
        
        # Run tleap command
        print("Running tleap to calculate system charge...")
        # Find tleap executable dynamically
        try:
            # First try to find tleap in PATH
            which_result = subprocess.run(['which', 'tleap'], capture_output=True, text=True)
            if which_result.returncode == 0:
                tleap_path = which_result.stdout.strip()
            else:
                # Fallback: try common conda environment paths
                conda_env = os.environ.get('CONDA_DEFAULT_ENV', 'MD_pipeline')
                conda_prefix = os.environ.get('CONDA_PREFIX', '')
                if conda_prefix:
                    tleap_path = os.path.join(conda_prefix, 'bin', 'tleap')
                else:
                    # Last resort: assume it's in PATH
                    tleap_path = 'tleap'
            
            cmd = f"{tleap_path} -f calc_charge_on_system.in"
            result = subprocess.run(cmd, shell=True, cwd=str(OUTPUT_DIR), capture_output=True, text=True)
        except Exception as e:
            # Fallback to simple tleap command
            cmd = f"tleap -f calc_charge_on_system.in"
            result = subprocess.run(cmd, shell=True, cwd=str(OUTPUT_DIR), capture_output=True, text=True)
        
        print(f"tleap return code: {result.returncode}")
        print(f"tleap stdout: {result.stdout}")
        print(f"tleap stderr: {result.stderr}")
        
        # Check if we got the charge information even if tleap had a non-zero exit code
        # (tleap often returns non-zero when run non-interactively but still calculates charge)
        if 'Total unperturbed charge' not in result.stdout and 'Total charge' not in result.stdout:
            return jsonify({'error': f'tleap failed to calculate charge. Error: {result.stderr}'}), 500
        
        # Parse the output to find the net charge
        output_lines = result.stdout.split('\n')
        net_charge = None
        
        for line in output_lines:
            if 'Total unperturbed charge' in line or 'Total charge' in line:
                # Look for patterns like "Total charge: -3.0000" or "Total unperturbed charge: -3.0000"
                import re
                charge_match = re.search(r'charge[:\s]+(-?\d+\.?\d*)', line)
                if charge_match:
                    net_charge = float(charge_match.group(1))
                    break
        
        if net_charge is None:
            return jsonify({'error': 'Could not extract net charge from tleap output'}), 500
        
        # Suggest ion addition
        if net_charge > 0:
            suggestion = f"Add {int(net_charge)} Cl- ions to neutralize the system"
            ion_type = "Cl-"
            ion_count = int(net_charge)
        elif net_charge < 0:
            suggestion = f"Add {int(abs(net_charge))} Na+ ions to neutralize the system"
            ion_type = "Na+"
            ion_count = int(abs(net_charge))
        else:
            suggestion = "System is already neutral, no ions needed"
            ion_type = "None"
            ion_count = 0
        
        return jsonify({
            'success': True,
            'net_charge': net_charge,
            'suggestion': suggestion,
            'ion_type': ion_type,
            'ion_count': ion_count,
            'ligand_present': ligand_present
        })
        
    except Exception as e:
        logger.error(f"Error calculating net charge: {str(e)}")
        return jsonify({'error': f'Internal server error: {str(e)}'}), 500

@app.route('/api/generate-all-files', methods=['POST'])
def generate_all_files():
    """Generate all simulation input files based on UI parameters"""
    try:
        data = request.get_json()
        
        # Get simulation parameters from UI
        cutoff_distance = data.get('cutoff_distance', 10.0)
        temperature = data.get('temperature', 310.0)
        pressure = data.get('pressure', 1.0)
        
        # Get step parameters
        restrained_steps = data.get('restrained_steps', 10000)
        restrained_force = data.get('restrained_force', 10.0)
        min_steps = data.get('min_steps', 20000)
        npt_heating_steps = data.get('npt_heating_steps', 50000)
        npt_equilibration_steps = data.get('npt_equilibration_steps', 100000)
        production_steps = data.get('production_steps', 1000000)
        # Integration time step (ps)
        dt = data.get('timestep', 0.002)
        
        # Get force field parameters
        force_field = data.get('force_field', 'ff14SB')
        water_model = data.get('water_model', 'TIP3P')
        add_ions = data.get('add_ions', 'None')
        distance = data.get('distance', 10.0)
        
        # Validation warnings
        warnings = []
        if restrained_steps < 5000:
            warnings.append("Restrained minimization steps should be at least 5000")
        if min_steps < 10000:
            warnings.append("Minimization steps should be at least 10000")
        
        # Count total residues in tleap_ready.pdb
        tleap_ready_file = OUTPUT_DIR / "tleap_ready.pdb"
        if not tleap_ready_file.exists():
            return jsonify({'error': 'tleap_ready.pdb not found. Please prepare structure first.'}), 400
        
        total_residues = count_residues_in_pdb(str(tleap_ready_file))
        
        # Generate min_restrained.in
        generate_min_restrained_file(restrained_steps, restrained_force, total_residues, cutoff_distance)
        
        # Generate min.in
        generate_min_file(min_steps, cutoff_distance)
        
        # Generate HeatNPT.in
        generate_heat_npt_file(npt_heating_steps, temperature, pressure, cutoff_distance, dt)
        
        # Generate mdin_equi.in (NPT Equilibration)
        generate_npt_equilibration_file(npt_equilibration_steps, temperature, pressure, cutoff_distance, dt)
        
        # Check if plumed.dat exists in output folder
        plumed_file = OUTPUT_DIR / 'plumed.dat'
        use_plumed = plumed_file.exists()
        
        # Generate mdin_prod.in (Production)
        generate_production_file(production_steps, temperature, pressure, cutoff_distance, dt, use_plumed=use_plumed)
        
        # Generate force field parameters
        ff_files_generated = []
        try:
            generate_ff_parameters_file(force_field, water_model, add_ions, distance)
            
            # Find tleap executable
            tleap_path = None
            try:
                result = subprocess.run(['which', 'tleap'], capture_output=True, text=True)
                if result.returncode == 0:
                    tleap_path = result.stdout.strip()
            except:
                pass
            
            if not tleap_path:
                conda_prefix = os.environ.get('CONDA_PREFIX')
                if conda_prefix:
                    tleap_path = os.path.join(conda_prefix, 'bin', 'tleap')
                else:
                    tleap_path = '/home/hn533621/.conda/envs/MD_pipeline/bin/tleap'
            
            # Run tleap to generate force field parameters
            cmd = f"{tleap_path} -f generate_ff_parameters.in"
            result = subprocess.run(cmd, shell=True, cwd=str(OUTPUT_DIR), 
                                  capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                warnings.append(f"Force field generation failed: {result.stderr}")
            else:
                # Check if key output files were created
                ff_output_files = ['protein.prmtop', 'protein.inpcrd', 'protein_solvated.pdb']
                for ff_file in ff_output_files:
                    if (OUTPUT_DIR / ff_file).exists():
                        ff_files_generated.append(ff_file)
                
                if len(ff_files_generated) == 0:
                    warnings.append("Force field parameter files were not generated")
                
        except Exception as ff_error:
            warnings.append(f"Force field generation error: {str(ff_error)}")
        
        # Generate PBS submit script into output
        pbs_generated = generate_submit_pbs_file(use_plumed=use_plumed)

        all_files = [
            'min_restrained.in',
            'min.in', 
            'HeatNPT.in',
            'mdin_equi.in',
            'mdin_prod.in'
        ] + ff_files_generated

        if pbs_generated:
            all_files.append('submit_job.pbs')
        
        return jsonify({
            'success': True,
            'message': f'All simulation files generated successfully ({len(all_files)} files)',
            'warnings': warnings,
            'files_generated': all_files
        })
        
    except Exception as e:
        logger.error(f"Error generating simulation files: {str(e)}")
        return jsonify({'error': f'Internal server error: {str(e)}'}), 500

def count_residues_in_pdb(pdb_file):
    """Count total number of residues in PDB file"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        residues = set()
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract residue number (columns 23-26)
                residue_num = line[22:26].strip()
                if residue_num:
                    residues.add(residue_num)
        
        return len(residues)
    except Exception as e:
        logger.error(f"Error counting residues: {str(e)}")
        return 607  # Default fallback

def generate_min_restrained_file(steps, force_constant, total_residues, cutoff):
    """Generate min_restrained.in file"""
    content = f"""initial minimization solvent + ions
 &cntrl
  imin   = 1,
  maxcyc = {steps},
  ncyc   = {steps // 2},
  ntb    = 1,
  ntr    = 1,
  ntxo   = 1,	
  cut    = {cutoff}
/
Restrain 
{force_constant}
RES 1 {total_residues}
END
END

"""
    
    with open(OUTPUT_DIR / "min_restrained.in", 'w') as f:
        f.write(content)

def generate_min_file(steps, cutoff):
    """Generate min.in file"""
    content = f"""Minimization
&cntrl
imin=1,
maxcyc={steps},
ncyc={steps // 4},
ntb=1,
cut={cutoff},
igb=0,
ntr=0,
/

"""
    
    with open(OUTPUT_DIR / "min.in", 'w') as f:
        f.write(content)

def generate_heat_npt_file(steps, temperature, pressure, cutoff, dt=0.002):
    """Generate HeatNPT.in file with temperature ramping"""
    # Calculate step divisions: 20%, 20%, 20%, 40%
    step1 = int(steps * 0.2)
    step2 = int(steps * 0.2)
    step3 = int(steps * 0.2)
    step4 = int(steps * 0.4)
    
    # Calculate temperature values: 3%, 66%, 100%
    temp1 = temperature * 0.03
    temp2 = temperature * 0.66
    temp3 = temperature
    temp4 = temperature
    
    content = f"""Heat
 &cntrl
  imin = 0, irest = 0, ntx = 1,
  ntb = 2, pres0 = {pressure}, ntp = 1,
  taup = 2.0,
  cut = {cutoff}, ntr = 0,
  ntc = 2, ntf = 2,
  tempi = 0, temp0 = {temperature},
  ntt = 3, gamma_ln = 1.0,
  nstlim = {steps}, dt = {dt},
  ntpr = 2000, ntwx = 2000, ntwr = 2000
 /
&wt type='TEMP0', istep1=0, istep2={step1}, value1=0.0, value2={temp1} /
&wt type='TEMP0', istep1={step1+1}, istep2={step1+step2}, value1={temp1}, value2={temp2} /
&wt type='TEMP0', istep1={step1+step2+1}, istep2={step1+step2+step3}, value1={temp2}, value2={temp3} /
&wt type='TEMP0', istep1={step1+step2+step3+1}, istep2={steps}, value1={temp3}, value2={temp4} /
&wt type='END' /

"""
    
    with open(OUTPUT_DIR / "HeatNPT.in", 'w') as f:
        f.write(content)

def generate_npt_equilibration_file(steps, temperature, pressure, cutoff, dt=0.002):
    """Generate mdin_equi.in file for NPT equilibration"""
    content = f"""NPT Equilibration
&cntrl
  imin=0,
  ntx=1,
  irest=0,
  pres0={pressure},
  taup=1.0,
  temp0={temperature},
  tempi={temperature},
  nstlim={steps},
  dt={dt},
  ntf=2,
  ntc=2,
  ntpr=500,
  ntwx=500,
  ntwr=500,
  cut={cutoff},
  ntb=2,
  ntp=1,
  ntt=3,
  gamma_ln=3.0,
  ig=-1,
  iwrap=1,
  ntr=0,
/

"""
    
    with open(OUTPUT_DIR / "mdin_equi.in", 'w') as f:
        f.write(content)

def generate_production_file(steps, temperature, pressure, cutoff, dt=0.002, use_plumed=False):
    """Generate mdin_prod.in file for production run"""
    content = f"""Production Run
&cntrl
  imin=0,
  ntx=1,
  irest=0,
  pres0={pressure},
  taup=1.0,
  temp0={temperature},
  tempi={temperature},
  nstlim={steps},
  dt={dt},
  ntf=2,
  ntc=2,
  ntpr=1000,
  ntwx=1000,
  ntwr=1000,
  cut={cutoff},
  ntb=2,
  ntp=1,
  ntt=3,
  gamma_ln=3.0,
  ig=-1,
  iwrap=1,
  ntr=0,
"""
    
    # Add PLUMED lines if plumed.dat exists
    if use_plumed:
        content += "  plumed=1,\n"
        content += "  plumedfile='plumed.dat'\n"
    
    content += "/\n\n"
    
    with open(OUTPUT_DIR / "mdin_prod.in", 'w') as f:
        f.write(content)

def generate_submit_pbs_file(use_plumed=False):
    """Generate submit_job.pbs file for SLURM job submission"""
    try:
        # Get absolute path to output directory
        output_dir_abs = OUTPUT_DIR.resolve()
        
        # Build PBS script content
        content = """#!/bin/bash
#SBATCH -D {working_dir}  # Critical: Sets working dir
#SBATCH --job-name=job_name
#SBATCH --partition=defq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=168:00:00


module load amber/24
""".format(working_dir=str(output_dir_abs))
        
        # Add PLUMED module if plumed.dat exists
        if use_plumed:
            content += "module load plumed/2.9.1\n"
        
        content += """
pmemd.cuda -O -i min_restrained.in -o min_restrained.out -p protein.prmtop -c protein.inpcrd -r min_res.ncrst -x min_res.nc -ref Leu15_Glu15.inpcrd -inf min_res.mdinfo
pmemd.cuda -O -i min.in -o min.out -p protein.prmtop -c min_res.ncrst -r min.ncrst -x min.nc -inf min.mdinfo
pmemd.cuda -O -i HeatNPT.in -o HeatNPT.out -p protein.prmtop -c min.ncrst -r HeatNPT.ncrst -x HeatNPT.nc -inf HeatNPT.mdinfo
pmemd.cuda -O -i mdin_equi.in -o mdin_equi.out -p protein.prmtop -c HeatNPT.ncrst -r mdin_equi.ncrst -x mdin_equi.nc -inf mdin_equi.mdinfo -ref HeatNPT.ncrst
pmemd.cuda -O -i mdin_prod.in -o mdin_prod.out -p protein.prmtop -c mdin_equi.ncrst -r mdin_prod.ncrst -x mdin_prod.nc -inf mdin_prod.mdinfo -ref mdin_equi.ncrst
"""
        
        # Write submit_job.pbs file
        with open(OUTPUT_DIR / "submit_job.pbs", 'w') as f:
            f.write(content)
        
        logger.info(f"Generated submit_job.pbs in {OUTPUT_DIR}")
        return True
    except Exception as e:
        logger.error(f"Error generating submit_job.pbs: {e}")
        return False

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({'status': 'healthy', 'message': 'MD Simulation Pipeline API is running'})

@app.route('/api/clean-output', methods=['POST'])
def clean_output():
    """Clean output folder endpoint"""
    try:
        print("DEBUG: clean-output endpoint called")
        if clean_and_create_output_folder():
            return jsonify({'success': True, 'message': 'Output folder cleaned successfully'})
        else:
            return jsonify({'success': False, 'error': 'Failed to clean output folder'}), 500
    except Exception as e:
        print(f"DEBUG: Error in clean-output: {str(e)}")
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/save-pdb-file', methods=['POST'])
def save_pdb_file():
    """Save PDB file to output directory"""
    try:
        data = request.get_json()
        pdb_content = data.get('pdb_content', '')
        filename = data.get('filename', 'input.pdb')
        
        if not pdb_content:
            return jsonify({'success': False, 'error': 'No PDB content provided'}), 400
        
        # Save to output directory as 0_original_input.pdb
        output_file = OUTPUT_DIR / "0_original_input.pdb"
        with open(output_file, 'w') as f:
            f.write(pdb_content)
        
        logger.info(f"Saved PDB file to {output_file}")
        return jsonify({
            'success': True,
            'message': f'PDB file saved successfully',
            'file_path': str(output_file)
        })
    except Exception as e:
        logger.error(f"Error saving PDB file: {str(e)}")
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/save-plumed-file', methods=['POST'])
def save_plumed_file():
    """Save PLUMED file to output directory"""
    try:
        data = request.get_json()
        plumed_content = data.get('plumed_content', '')
        filename = data.get('filename', 'plumed.dat')
        
        if not plumed_content:
            return jsonify({'success': False, 'error': 'No PLUMED content provided'}), 400
        
        # Ensure filename has .dat extension if not provided
        if not filename.endswith('.dat'):
            filename = filename if '.' in filename else f"{filename}.dat"
        
        # Save to output directory
        output_file = OUTPUT_DIR / filename
        with open(output_file, 'w') as f:
            f.write(plumed_content)
        
        logger.info(f"Saved PLUMED file to {output_file}")
        return jsonify({
            'success': True,
            'message': f'PLUMED file saved successfully to output/{filename}',
            'file_path': str(output_file),
            'filename': filename
        })
    except Exception as e:
        logger.error(f"Error saving PLUMED file: {str(e)}")
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/download-output-zip', methods=['GET'])
def download_output_zip():
    """Create a ZIP of the output folder and return it for download"""
    try:
        if not OUTPUT_DIR.exists():
            return jsonify({'error': 'Output directory not found'}), 404

        import tempfile
        import shutil

        # Create a temporary zip file
        tmp_dir = tempfile.mkdtemp()
        zip_base = os.path.join(tmp_dir, 'output')
        zip_path = shutil.make_archive(zip_base, 'zip', root_dir=str(OUTPUT_DIR))

        # Send file for download
        return send_file(zip_path, as_attachment=True, download_name='output.zip')
    except Exception as e:
        logger.error(f"Error creating output ZIP: {str(e)}")
        return jsonify({'error': f'Failed to create ZIP: {str(e)}'}), 500

@app.route('/api/get-generated-files', methods=['GET'])
def get_generated_files():
    """Return contents of known generated input files for preview"""
    try:
        files_to_read = [
            'min_restrained.in',
            'min.in',
            'HeatNPT.in',
            'mdin_equi.in',
            'mdin_prod.in',
            'submit_job.pbs'
        ]
        # Files to exclude from preview (intermediate/utility files)
        excluded_files = [
            'calc_charge_on_system.in',
            'generate_ff_parameters.in',
            'sqm.in'
        ]
        # Note: Force field parameter files (protein.prmtop, protein.inpcrd, protein_solvated.pdb) 
        # are excluded from preview as they are binary/large files
        
        # Also include any user-created .in files in the output directory
        user_created_files = []
        try:
            for file_path in OUTPUT_DIR.glob("*.in"):
                filename = file_path.name
                # Exclude standard files and utility files
                if filename not in files_to_read and filename not in excluded_files:
                    user_created_files.append(filename)
        except Exception as e:
            logger.warning(f"Error scanning for user-created files: {e}")
        
        # Combine standard files and user-created files
        all_files = files_to_read + sorted(user_created_files)
        
        result = {}
        for name in all_files:
            path = OUTPUT_DIR / name
            if path.exists():
                try:
                    with open(path, 'r') as f:
                        result[name] = f.read()
                except Exception as fe:
                    result[name] = f"<error reading file: {fe}>"
            else:
                result[name] = "<file not found>"
        return jsonify({'success': True, 'files': result})
    except Exception as e:
        logger.error(f"Error reading generated files: {str(e)}")
        return jsonify({'error': f'Failed to read files: {str(e)}'}), 500

@app.route('/api/save-file', methods=['POST'])
def save_file():
    """Save edited file content back to the output directory"""
    try:
        data = request.get_json()
        filename = data.get('filename')
        content = data.get('content')
        
        if not filename:
            return jsonify({'success': False, 'error': 'Filename is required'}), 400
        
        if content is None:
            return jsonify({'success': False, 'error': 'Content is required'}), 400
        
        # Security: Only allow saving files that are in the allowed list
        allowed_files = [
            'min_restrained.in',
            'min.in',
            'HeatNPT.in',
            'mdin_equi.in',
            'mdin_prod.in',
            'submit_job.pbs'
        ]
        
        if filename not in allowed_files:
            return jsonify({'success': False, 'error': f'File "{filename}" is not allowed to be edited'}), 403
        
        # Prevent directory traversal attacks
        if '/' in filename or '\\' in filename or '..' in filename:
            return jsonify({'success': False, 'error': 'Invalid filename'}), 400
        
        # Write file
        file_path = OUTPUT_DIR / filename
        try:
            with open(file_path, 'w') as f:
                f.write(content)
            
            logger.info(f"File {filename} saved successfully")
            return jsonify({'success': True, 'message': f'File {filename} saved successfully'})
        except Exception as e:
            logger.error(f"Error writing file {filename}: {str(e)}")
            return jsonify({'success': False, 'error': f'Failed to write file: {str(e)}'}), 500
            
    except Exception as e:
        logger.error(f"Error saving file: {str(e)}")
        return jsonify({'success': False, 'error': f'Internal server error: {str(e)}'}), 500

@app.route('/api/save-new-file', methods=['POST'])
def save_new_file():
    """Save a new simulation file created by the user"""
    try:
        data = request.get_json()
        filename = data.get('filename')
        content = data.get('content')
        
        if not filename:
            return jsonify({'success': False, 'error': 'Filename is required'}), 400
        
        if content is None:
            return jsonify({'success': False, 'error': 'Content is required'}), 400
        
        # Validate filename - must end with .in
        if not filename.endswith('.in'):
            return jsonify({'success': False, 'error': 'File name must end with .in extension'}), 400
        
        # Prevent directory traversal attacks
        if '/' in filename or '\\' in filename or '..' in filename:
            return jsonify({'success': False, 'error': 'Invalid filename'}), 400
        
        # Write file
        file_path = OUTPUT_DIR / filename
        try:
            with open(file_path, 'w') as f:
                f.write(content)
            
            logger.info(f"New file {filename} saved successfully")
            return jsonify({'success': True, 'message': f'File {filename} saved successfully'})
        except Exception as e:
            logger.error(f"Error writing new file {filename}: {str(e)}")
            return jsonify({'success': False, 'error': f'Failed to write file: {str(e)}'}), 500
            
    except Exception as e:
        logger.error(f"Error saving new file: {str(e)}")
        return jsonify({'success': False, 'error': f'Internal server error: {str(e)}'}), 500

def get_ligand_residue_name():
    """Extract first ligand residue name from tleap_ready.pdb (for backward compatibility)"""
    ligand_names = get_all_ligand_residue_names()
    return ligand_names[0] if ligand_names else "LIG"

def generate_ff_parameters_file(force_field, water_model, add_ions, distance):
    """Generate the final force field parameters file with dynamic values"""
    # Debug logging
    print(f"DEBUG: force_field={force_field}, water_model={water_model}, add_ions={add_ions}, distance={distance}")
    
    # Check if ligands are present - look for residue-named files first, then fallback to numbered files
    ligand_mol2_files = []
    ligand_frcmod_files = []
    ligand_resname_map = {}  # Maps residue name to (mol2_file, frcmod_file)
    
    # First, try to find residue-named files (e.g., O9C.mol2, O9C.frcmod)
    unique_resnames = get_all_ligand_residue_names()
    for resname in unique_resnames:
        mol2_file = OUTPUT_DIR / f"{resname}.mol2"
        frcmod_file = OUTPUT_DIR / f"{resname}.frcmod"
        if mol2_file.exists() and frcmod_file.exists():
            ligand_resname_map[resname] = (mol2_file, frcmod_file)
            ligand_mol2_files.append(mol2_file)
            ligand_frcmod_files.append(frcmod_file)
    
    # Fallback: check for numbered files (backward compatibility)
    if not ligand_mol2_files:
        numbered_mol2 = sorted(OUTPUT_DIR.glob("4_ligands_corrected_*.mol2"))
        numbered_frcmod = sorted(OUTPUT_DIR.glob("4_ligands_corrected_*.frcmod"))
        if numbered_mol2 and numbered_frcmod:
            ligand_mol2_files = numbered_mol2
            ligand_frcmod_files = numbered_frcmod
            # Try to map to residue names
            resnames = get_all_ligand_residue_names()
            for i, (mol2_file, frcmod_file) in enumerate(zip(ligand_mol2_files, ligand_frcmod_files)):
                # Extract residue name from mol2 file if possible
                resname = get_residue_name_from_mol2(mol2_file) if mol2_file.exists() else None
                if not resname:
                    # Try to get from tleap_ready.pdb
                    if resnames and i < len(resnames):
                        resname = resnames[i]
                    else:
                        resname = f"LIG{len(ligand_resname_map) + 1}"
                # Only add if not already in map (avoid duplicates)
                if resname not in ligand_resname_map:
                    ligand_resname_map[resname] = (mol2_file, frcmod_file)
    
    # Final fallback: single ligand file (backward compatibility)
    if not ligand_mol2_files:
        single_mol2 = OUTPUT_DIR / "4_ligands_corrected.mol2"
        single_frcmod = OUTPUT_DIR / "4_ligands_corrected.frcmod"
        if single_mol2.exists() and single_frcmod.exists():
            ligand_mol2_files = [single_mol2]
            ligand_frcmod_files = [single_frcmod]
            resnames = get_all_ligand_residue_names()
            if resnames:
                ligand_resname_map[resnames[0]] = (single_mol2, single_frcmod)
            else:
                ligand_resname_map["LIG"] = (single_mol2, single_frcmod)
    
    ligand_present = len(ligand_mol2_files) > 0 and len(ligand_frcmod_files) > 0
    
    # Build the content dynamically
    content = f"source leaprc.protein.{force_field}\n"
    
    # Add water model source
    print(f"DEBUG: water_model={water_model}")
    if water_model.lower() == "tip3p":
        content += "source leaprc.water.tip3p\n"
    elif water_model == "spce":
        content += "source leaprc.water.spce\n"
    
    # Add ligand-related commands only if ligands are present
    if ligand_present:
        content += "source leaprc.gaff2\n\n"
        
        # Load each unique ligand parameter and structure only once
        # Use sorted to ensure consistent ordering
        for resname in sorted(ligand_resname_map.keys()):
            mol2_file, frcmod_file = ligand_resname_map[resname]
            content += f"loadamberparams {frcmod_file.name}\n"
            content += f"{resname} = loadmol2 {mol2_file.name}\n"
        content += "\n"
    else:
        content += "\n"
    
    content += "x = loadpdb tleap_ready.pdb\n\n"
    content += "charge x\n\n"
    
    # Add ions based on selection
    if add_ions == "Na+":
        content += "addions x Na+ 0.0\n\n"
    elif add_ions == "Cl-":
        content += "addions x Cl- 0.0\n\n"
    # If "None", skip adding ions
    
    # Add solvation with selected water model and distance
    if water_model.lower() == "tip3p":
        content += f"solvateBox x TIP3PBOX {distance}\n\n"
    elif water_model.lower() == "spce":
        content += f"solvateBox x SPCBOX {distance}\n\n"
    
    content += "saveamberparm x protein.prmtop protein.inpcrd\n\n"
    content += "savepdb x protein_solvated.pdb\n\n"
    content += "quit\n"
    
    # Debug: print the generated content
    print("DEBUG: Generated content:")
    print(content)
    
    # Write the file
    with open(OUTPUT_DIR / "generate_ff_parameters.in", 'w') as f:
        f.write(content)

def get_residue_name_from_pdb(pdb_file):
    """Extract residue name from a ligand PDB file"""
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    # Extract residue name (columns 18-20)
                    residue_name = line[17:20].strip()
                    if residue_name and residue_name not in ['HOH', 'WAT', 'TIP', 'SPC', 'NA', 'CL']:
                        return residue_name
        return None
    except Exception as e:
        logger.warning(f"Could not extract residue name from {pdb_file}: {e}")
        return None

def get_residue_name_from_mol2(mol2_file):
    """Extract residue name from a mol2 file (from @<TRIPOS>MOLECULE section)"""
    try:
        with open(mol2_file, 'r') as f:
            lines = f.readlines()
            # Find @<TRIPOS>MOLECULE section
            in_molecule = False
            for i, line in enumerate(lines):
                if '@<TRIPOS>MOLECULE' in line:
                    in_molecule = True
                    # The next line is the molecule name/residue name
                    if i + 1 < len(lines):
                        resname = lines[i + 1].strip()
                        # Remove any extra whitespace or comments
                        resname = resname.split()[0] if resname.split() else resname
                        return resname
        return None
    except Exception as e:
        logger.warning(f"Could not extract residue name from {mol2_file}: {e}")
        return None

def get_all_ligand_residue_names():
    """Extract all unique ligand residue names from tleap_ready.pdb"""
    ligand_names = []
    try:
        tleap_ready_path = OUTPUT_DIR / "tleap_ready.pdb"
        if not tleap_ready_path.exists():
            return []
        
        seen_residues = set()
        with open(tleap_ready_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    # Extract residue name (columns 18-20)
                    residue_name = line[17:20].strip()
                    if residue_name and residue_name not in ['HOH', 'WAT', 'TIP', 'SPC', 'NA', 'CL']:
                        if residue_name not in seen_residues:
                            ligand_names.append(residue_name)
                            seen_residues.add(residue_name)
        
        return ligand_names
    except Exception as e:
        logger.warning(f"Could not extract ligand residue names: {e}")
        return []

@app.route('/api/generate-ff-parameters', methods=['POST'])
def generate_ff_parameters():
    """Generate final force field parameters using tleap"""
    try:
        data = request.get_json()
        force_field = data.get('force_field', 'ff14SB')
        water_model = data.get('water_model', 'TIP3P')
        add_ions = data.get('add_ions', 'None')
        distance = data.get('distance', 10.0)
        
        # Generate the dynamic input file
        generate_ff_parameters_file(force_field, water_model, add_ions, distance)
        
        # Find tleap executable
        tleap_path = None
        try:
            result = subprocess.run(['which', 'tleap'], capture_output=True, text=True)
            if result.returncode == 0:
                tleap_path = result.stdout.strip()
        except:
            pass
        
        if not tleap_path:
            conda_prefix = os.environ.get('CONDA_PREFIX')
            if conda_prefix:
                tleap_path = os.path.join(conda_prefix, 'bin', 'tleap')
            else:
                tleap_path = '/home/hn533621/.conda/envs/MD_pipeline/bin/tleap'
        
        # Run tleap
        cmd = f"{tleap_path} -f generate_ff_parameters.in"
        result = subprocess.run(cmd, shell=True, cwd=str(OUTPUT_DIR), 
                              capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            logger.error(f"tleap failed: {result.stderr}")
            return jsonify({
                'success': False, 
                'error': f'tleap failed: {result.stderr}'
            }), 500
        
        # Check if key output files were created
        output_files = ['protein.prmtop', 'protein.inpcrd', 'protein_solvated.pdb']
        missing_files = [f for f in output_files if not (OUTPUT_DIR / f).exists()]
        
        if missing_files:
            return jsonify({
                'success': False,
                'error': f'Missing output files: {", ".join(missing_files)}'
            }), 500
        
        return jsonify({
            'success': True,
            'message': 'Force field parameters generated successfully',
            'files_generated': output_files
        })
        
    except subprocess.TimeoutExpired:
        return jsonify({
            'success': False,
            'error': 'tleap command timed out after 5 minutes'
        }), 500
    except Exception as e:
        logger.error(f"Error generating FF parameters: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Failed to generate force field parameters: {str(e)}'
        }), 500

@app.route('/api/detect-missing-residues', methods=['POST'])
def detect_missing_residues_endpoint():
    """Detect missing residues in the loaded PDB structure"""
    try:
        # Check if original input file exists
        original_pdb_path = OUTPUT_DIR / "0_original_input.pdb"
        if not original_pdb_path.exists():
            return jsonify({
                'success': False,
                'error': 'No PDB file loaded. Please load a PDB file first.'
            }), 400
        
        # Get PDB ID from the file
        try:
            pdb_id = get_pdb_id_from_pdb_file(str(original_pdb_path))
        except ValueError as e:
            return jsonify({
                'success': False,
                'error': f'Could not determine PDB ID: {str(e)}'
            }), 400
        
        # Detect missing residues
        missing = detect_missing_residues(pdb_id)
        
        # Get chain sequences
        chain_sequences = get_chain_sequences(pdb_id)
        
        # Find chains with missing residues that have sequences available
        chains_with_missing = {
            chain: chain_sequences[chain]
            for chain in missing
            if chain in chain_sequences
        }
        
        # Format missing residues info for display
        missing_info = {}
        for chain, missing_list in missing.items():
            missing_info[chain] = {
                'count': len(missing_list),
                'residues': missing_list
            }
        
        # Get first residue number for each chain from the PDB file
        # Also calculate the starting residue number for the sequence viewer
        # (accounting for missing residues before the first PDB residue)
        chain_first_residue = {}
        chain_sequence_start = {}
        try:
            original_pdb_path = OUTPUT_DIR / "0_original_input.pdb"
            if original_pdb_path.exists():
                with open(original_pdb_path, 'r') as f:
                    pdb_lines = f.readlines()
                    
                # First pass: find first residue number for each chain
                for line in pdb_lines:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        chain_id = line[21:22].strip()
                        if chain_id and chain_id not in chain_first_residue:
                            # Extract residue number (columns 22-26, but we need to handle insertion codes)
                            residue_str = line[22:26].strip()
                            try:
                                # Try to extract just the number part (handle negative numbers)
                                import re
                                match = re.match(r'(-?\d+)', residue_str)
                                if match:
                                    residue_num = int(match.group(1))
                                    chain_first_residue[chain_id] = residue_num
                            except:
                                pass
                
                # Second pass: calculate sequence start for each chain
                # We want to find the first residue number that should be displayed
                # This is the first PDB residue minus the count of missing residues before it
                # Example: If PDB starts at 189 and residues 173-188 are missing (16 residues),
                # then sequence_start = 189 - 16 = 173
                for chain_id, first_pdb_residue in chain_first_residue.items():
                    # Find the minimum missing residue number before first_pdb_residue
                    # This tells us where the sequence should start displaying
                    min_missing_before = None
                    if chain_id in missing_info:
                        for resname, resnum in missing_info[chain_id]['residues']:
                            if resnum < first_pdb_residue:
                                if min_missing_before is None or resnum < min_missing_before:
                                    min_missing_before = resnum
                    
                    if min_missing_before is not None:
                        # Sequence should start from the first missing residue before PDB start
                        # This accounts for all missing residues before the first PDB residue
                        sequence_start = min_missing_before
                    else:
                        # No missing residues before first PDB residue, start from first PDB residue
                        sequence_start = first_pdb_residue
                    
                    chain_sequence_start[chain_id] = sequence_start
        except Exception as e:
            logger.warning(f"Could not determine first residue numbers: {str(e)}")
        
        return jsonify({
            'success': True,
            'pdb_id': pdb_id,
            'missing_residues': missing_info,
            'chains_with_missing': list(chains_with_missing.keys()),
            'chain_sequences': chain_sequences,
            'chain_first_residue': chain_first_residue,
            'chain_sequence_start': chain_sequence_start
        })
        
    except Exception as e:
        logger.error(f"Error detecting missing residues: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Failed to detect missing residues: {str(e)}'
        }), 500

@app.route('/api/trim-residues', methods=['POST'])
def trim_residues_endpoint():
    """Trim residues from edges of chain sequences"""
    try:
        data = request.get_json()
        chain_sequences = data.get('chain_sequences', {})
        trim_specs = data.get('trim_specs', {})
        pdb_id = data.get('pdb_id')
        
        if not chain_sequences:
            return jsonify({
                'success': False,
                'error': 'No chain sequences provided'
            }), 400
        
        if not trim_specs:
            return jsonify({
                'success': False,
                'error': 'No trim specifications provided'
            }), 400
        
        # Apply trimming
        try:
            trimmed_sequences = trim_chains_sequences(chain_sequences, trim_specs)
        except ValueError as e:
            return jsonify({
                'success': False,
                'error': str(e)
            }), 400
        
        # Optionally write trimmed FASTA file if pdb_id is provided
        if pdb_id:
            try:
                write_fasta_for_missing_chains(
                    pdb_id, 
                    trimmed_sequences, 
                    output_dir=str(OUTPUT_DIR)
                )
                logger.info(f"Wrote trimmed FASTA file for PDB {pdb_id}")
            except Exception as e:
                logger.warning(f"Could not write trimmed FASTA file: {str(e)}")
        
        # Calculate trim info for response
        trim_info = {}
        for chain, spec in trim_specs.items():
            original_len = len(chain_sequences.get(chain, ''))
            trimmed_len = len(trimmed_sequences.get(chain, ''))
            trim_info[chain] = {
                'original_length': original_len,
                'trimmed_length': trimmed_len,
                'n_terminal_trimmed': spec.get('n_terminal', 0),
                'c_terminal_trimmed': spec.get('c_terminal', 0)
            }
        
        return jsonify({
            'success': True,
            'trimmed_sequences': trimmed_sequences,
            'trim_info': trim_info,
            'message': f'Successfully trimmed residues from {len(trim_specs)} chain(s)'
        })
        
    except Exception as e:
        logger.error(f"Error trimming residues: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Failed to trim residues: {str(e)}'
        }), 500

@app.route('/api/build-completed-structure', methods=['POST'])
def build_completed_structure_endpoint():
    """Build completed structure using ESMFold for selected chains"""
    try:
        data = request.get_json()
        selected_chains = data.get('selected_chains', [])
        
        if not selected_chains:
            return jsonify({
                'success': False,
                'error': 'No chains selected for completion'
            }), 400
        
        # Check if original input file exists
        original_pdb_path = OUTPUT_DIR / "0_original_input.pdb"
        if not original_pdb_path.exists():
            return jsonify({
                'success': False,
                'error': 'No PDB file loaded. Please load a PDB file first.'
            }), 400
        
        # Get PDB ID
        try:
            pdb_id = get_pdb_id_from_pdb_file(str(original_pdb_path))
        except ValueError as e:
            return jsonify({
                'success': False,
                'error': f'Could not determine PDB ID: {str(e)}'
            }), 400
        
        # Get chain sequences (use provided sequences if available, otherwise fetch)
        provided_sequences = data.get('chain_sequences', None)
        if provided_sequences:
            chain_sequences = provided_sequences
            logger.info("Using provided chain sequences (may be trimmed)")
        else:
            chain_sequences = get_chain_sequences(pdb_id)
        
        # Verify selected chains have sequences
        chains_to_process = []
        for chain in selected_chains:
            if chain in chain_sequences:
                chains_to_process.append(chain)
            else:
                logger.warning(f"Chain {chain} not found in chain sequences")
        
        if not chains_to_process:
            return jsonify({
                'success': False,
                'error': 'None of the selected chains have sequences available'
            }), 400
        
        # Create dictionary of chains with their sequences for FASTA writing
        chains_with_missing = {
            chain: chain_sequences[chain]
            for chain in chains_to_process
        }
        
        # Write FASTA file for the selected chains
        try:
            write_fasta_for_missing_chains(pdb_id, chains_with_missing, output_dir=str(OUTPUT_DIR))
            logger.info(f"Wrote FASTA file for chains: {chains_to_process}")
        except Exception as e:
            logger.warning(f"Could not write FASTA file: {str(e)}")
            # Don't fail the entire operation if FASTA writing fails
        
        # Run ESMFold for each selected chain
        esmfold_results = {}
        for chain in chains_to_process:
            logger.info(f"Running ESMFold for chain {chain}")
            seq = chain_sequences[chain]
            try:
                pdb_text = run_esmfold(seq)
                esmfold_results[chain] = pdb_text
                
                # Save each chain's ESMFold result
                esm_pdb_filename = OUTPUT_DIR / f"{pdb_id}_chain_{chain}_esmfold.pdb"
                with open(esm_pdb_filename, 'w') as f:
                    f.write(pdb_text)
                logger.info(f"Saved ESMFold result for chain {chain} to {esm_pdb_filename}")
            except Exception as e:
                logger.error(f"Error running ESMFold for chain {chain}: {str(e)}")
                return jsonify({
                    'success': False,
                    'error': f'ESMFold failed for chain {chain}: {str(e)}'
                }), 500
        
        # Rebuild PDB using PyMOL
        output_pdb = OUTPUT_DIR / "0_complete_structure.pdb"
        try:
            # Use subprocess to run PyMOL in a separate process to avoid conflicts
            import tempfile
            import os
            
            # Create a standalone script that runs PyMOL operations
            script_content = f"""#!/usr/bin/env python3
import sys
import os

# Add parent directory to path
sys.path.insert(0, r'{str(Path(__file__).parent.parent)}')

# Change to output directory
os.chdir(r'{str(OUTPUT_DIR)}')

# Import and run rebuild
from Fill_missing_residues import rebuild_pdb_with_esmfold

try:
    rebuild_pdb_with_esmfold(
        r'{pdb_id}',
        {repr(chains_to_process)},
        output_pdb=r'{output_pdb.name}',
        original_pdb_path=r'{Path(original_pdb_path).name}'
    )
    print("SUCCESS: Rebuild completed")
except Exception as e:
    print(f"ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
"""
            
            # Write script to temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as script_file:
                script_file.write(script_content)
                script_path = script_file.name
            
            try:
                # Make script executable
                os.chmod(script_path, 0o755)
                
                # Run script in subprocess
                result = subprocess.run(
                    [sys.executable, script_path],
                    capture_output=True,
                    text=True,
                    timeout=300,
                    cwd=str(OUTPUT_DIR)
                )
                
                if result.returncode != 0:
                    error_msg = result.stderr or result.stdout
                    logger.error(f"PyMOL rebuild failed: {error_msg}")
                    # Check if it's a PyMOL initialization issue
                    if "pymol" in error_msg.lower() or "import" in error_msg.lower():
                        raise Exception(f"PyMOL initialization failed. Make sure PyMOL is installed and accessible. Error: {error_msg}")
                    raise Exception(f"Rebuild failed: {error_msg}")
                
                if "ERROR:" in result.stdout:
                    error_line = [line for line in result.stdout.split('\\n') if 'ERROR:' in line]
                    if error_line:
                        raise Exception(error_line[0].replace('ERROR:', '').strip())
                
                if not output_pdb.exists():
                    raise Exception("Output file was not created")
                
                logger.info(f"Completed structure saved to {output_pdb}")
                
            finally:
                # Clean up temporary script
                try:
                    os.unlink(script_path)
                except:
                    pass
                    
        except subprocess.TimeoutExpired:
            logger.error("PyMOL rebuild timed out after 5 minutes")
            return jsonify({
                'success': False,
                'error': 'PyMOL rebuild timed out. The structure might be too large. Please try again.'
            }), 500
        except Exception as e:
            logger.error(f"Error rebuilding PDB: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return jsonify({
                'success': False,
                'error': f'Failed to rebuild structure: {str(e)}'
            }), 500
        
        # Read the completed structure
        with open(output_pdb, 'r') as f:
            completed_content = f.read()
        
        return jsonify({
            'success': True,
            'message': f'Successfully completed structure for chains: {", ".join(chains_to_process)}',
            'completed_chains': chains_to_process,
            'completed_structure': completed_content
        })
        
    except Exception as e:
        logger.error(f"Error building completed structure: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Failed to build completed structure: {str(e)}'
        }), 500

@app.route('/api/get-completed-structure', methods=['GET'])
def get_completed_structure():
    """Get the completed structure PDB file if it exists"""
    try:
        completed_pdb_path = OUTPUT_DIR / "0_complete_structure.pdb"
        if not completed_pdb_path.exists():
            return jsonify({
                'success': False,
                'exists': False,
                'error': 'Completed structure not found'
            }), 404
        
        with open(completed_pdb_path, 'r') as f:
            content = f.read()
        
        return jsonify({
            'success': True,
            'exists': True,
            'content': content
        })
    except Exception as e:
        logger.error(f"Error reading completed structure: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/get-file', methods=['GET'])
def get_file():
    """Get a file from the output directory"""
    try:
        filename = request.args.get('filename')
        if not filename:
            return jsonify({
                'success': False,
                'error': 'Filename parameter required'
            }), 400
        
        # Security: only allow files from output directory
        file_path = OUTPUT_DIR / filename
        
        # Prevent directory traversal
        if not str(file_path).startswith(str(OUTPUT_DIR)):
            return jsonify({
                'success': False,
                'error': 'Invalid file path'
            }), 400
        
        if not file_path.exists():
            return jsonify({
                'success': False,
                'error': f'File {filename} not found'
            }), 404
        
        # Read file content
        with open(file_path, 'r') as f:
            content = f.read()
        
        return content, 200, {'Content-Type': 'text/plain'}
    except Exception as e:
        logger.error(f"Error reading file {filename}: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

if __name__ == '__main__':
    print("🧬 MD Simulation Pipeline")
    print("=========================")
    print("🌐 Starting Flask server...")
    print("📡 Backend API: http://localhost:5000")
    print("🔗 Web Interface: http://localhost:5000")
    print("")
    print("Press Ctrl+C to stop the server")
    print("")
    
    # Clean and create fresh output folder on startup
    print("🧹 Cleaning output folder...")
    clean_and_create_output_folder()
    print("✅ Output folder ready!")
    print("")
    
    app.run(debug=False, host='0.0.0.0', port=5000)
