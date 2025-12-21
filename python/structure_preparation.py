#!/usr/bin/env python3
"""
AMBER Structure Preparation Script using MDAnalysis
Complete pipeline: extract protein, add caps, handle ligands
"""

import os
import subprocess
import sys
import shutil
import logging

logger = logging.getLogger(__name__)

def run_command(cmd, description=""):
    """Run a command and return success status"""
    try:
        print(f"Running: {description}")
        print(f"Command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        print(f"Return code: {result.returncode}")
        if result.stdout:
            print(f"STDOUT: {result.stdout}")
        if result.stderr:
            print(f"STDERR: {result.stderr}")
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"Timeout: {description}")
        return False
    except Exception as e:
        print(f"Error running {description}: {str(e)}")
        return False

def extract_protein_only(pdb_content, output_file, selected_chains=None):
    """Extract protein without hydrogens using MDAnalysis. Optionally restrict to selected chains."""
    # Write input content to output file first
    with open(output_file, 'w') as f:
        f.write(pdb_content)
    
    try:
        # Run MDAnalysis command with the output file as input
        chain_sel = ''
        if selected_chains:
            chain_filters = ' or '.join([f'chain {c}' for c in selected_chains])
            chain_sel = f' and ({chain_filters})'
        selection = f"protein{chain_sel} and not name H* 1H* 2H* 3H*"
        cmd = f'python -c "import MDAnalysis as mda; u=mda.Universe(\'{output_file}\'); u.select_atoms(\'{selection}\').write(\'{output_file}\')"'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            raise Exception(f"MDAnalysis error: {result.stderr}")
        
        return True
    except Exception as e:
        print(f"Error in extract_protein_only: {e}")
        return False

def add_capping_groups(input_file, output_file):
    """Add ACE and NME capping groups using add_caps.py"""
    # First add caps
    temp_capped = output_file.replace('.pdb', '_temp.pdb')
    cmd = f"python add_caps.py -i {input_file} -o {temp_capped}"
    if not run_command(cmd, f"Adding capping groups to {input_file}"):
        return False
    
    # Then add TER cards using awk
    cmd = f"awk '/NME/{{nme=NR}} /ACE/ && nme && NR > nme {{print \"TER\"; nme=0}} {{print}}' {temp_capped} > {output_file}"
    if not run_command(cmd, f"Adding TER cards to {temp_capped}"):
        return False
    
    # Clean up temp file
    if os.path.exists(temp_capped):
        os.remove(temp_capped)
    
    return True

def extract_selected_chains(pdb_content, output_file, selected_chains):
    """Extract selected chains using PyMOL commands"""
    try:
        # Write input content to temp file
        temp_input = output_file.replace('.pdb', '_temp_input.pdb')
        with open(temp_input, 'w') as f:
            f.write(pdb_content)
        
        # Build chain selection string
        chain_filters = ' or '.join([f'chain {c}' for c in selected_chains])
        selection = f"({chain_filters}) and polymer.protein"
        
        # Use PyMOL to extract chains
        cmd = f'''python -c "
import pymol
pymol.finish_launching(['pymol', '-c'])
pymol.cmd.load('{temp_input}')
pymol.cmd.save('{output_file}', '{selection}')
pymol.cmd.quit()
"'''
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        
        # Clean up temp file
        if os.path.exists(temp_input):
            os.remove(temp_input)
        
        if result.returncode != 0:
            print(f"PyMOL chain extraction error: {result.stderr}")
            return False
        
        return True
    except Exception as e:
        print(f"Error extracting selected chains: {e}")
        return False

def extract_selected_ligands(pdb_content, output_file, selected_ligands):
    """Extract selected ligands using PyMOL commands"""
    try:
        # Write input content to temp file
        temp_input = output_file.replace('.pdb', '_temp_input.pdb')
        with open(temp_input, 'w') as f:
            f.write(pdb_content)
        
        # Build ligand selection string
        parts = []
        for lig in selected_ligands:
            resn = lig.get('resn', '').strip()
            chain = lig.get('chain', '').strip()
            if resn and chain:
                parts.append(f"(resn {resn} and chain {chain})")
            elif resn:
                parts.append(f"resn {resn}")
        
        if not parts:
            # No ligands to extract
            with open(output_file, 'w') as f:
                f.write('\n')
            return True
        
        selection = ' or '.join(parts)
        
        # Use PyMOL to extract ligands
        cmd = f'''python -c "
import pymol
pymol.finish_launching(['pymol', '-c'])
pymol.cmd.load('{temp_input}')
pymol.cmd.save('{output_file}', '{selection}')
pymol.cmd.quit()
"'''
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        
        # Clean up temp file
        if os.path.exists(temp_input):
            os.remove(temp_input)
        
        if result.returncode != 0:
            print(f"PyMOL ligand extraction error: {result.stderr}")
            return False
        
        return True
    except Exception as e:
        print(f"Error extracting selected ligands: {e}")
        return False

def extract_ligands(pdb_content, output_file, ligand_residue_name=None, selected_ligands=None):
    """Extract ligands using MDAnalysis. Optionally restrict to selected ligands (list of dicts with resn, chain, resi)."""
    # Write input content to output file first
    with open(output_file, 'w') as f:
        f.write(pdb_content)
    
    try:
        # Run MDAnalysis command with the output file as input
        if selected_ligands:
            # Build selection from provided ligand list (RESN-CHAIN groups)
            parts = []
            for lig in selected_ligands:
                resn = lig.get('resn', '').strip()
                chain = lig.get('chain', '').strip()
                if resn and chain:
                    parts.append(f"(resname {resn} and segid {chain})")
                elif resn:
                    parts.append(f"resname {resn}")
            if parts:
                selection = ' or '.join(parts)
                cmd = f'''python -c "
import MDAnalysis as mda
u = mda.Universe('{output_file}')
u.select_atoms('{selection}').write('{output_file}')
"'''
            else:
                cmd = f"python -c \"open('{output_file}','w').write('\\n')\""
        elif ligand_residue_name:
            # Use specified ligand residue name - extract from both ATOM and HETATM records
            cmd = f'''python -c "
import MDAnalysis as mda
u = mda.Universe('{output_file}')
# Extract specific ligand residue from both ATOM and HETATM records
u.select_atoms('resname {ligand_residue_name}').write('{output_file}')
"'''
        else:
            # Auto-detect ligand residues
            cmd = f'''python -c "
import MDAnalysis as mda
u = mda.Universe('{output_file}')
# Get all unique residue names from HETATM records
hetatm_residues = set()
for atom in u.atoms:
    if atom.record_type == 'HETATM':
        hetatm_residues.add(atom.resname)
# Remove water and ions
ligand_residues = hetatm_residues - {{'HOH', 'WAT', 'TIP3', 'TIP4', 'SPC', 'SPCE', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 'CD', 'HG', 'PB', 'SR', 'BA', 'RB', 'CS', 'LI', 'F', 'BR', 'I', 'PO4', 'PO3', 'H2PO4', 'HPO4', 'H3PO4', 'SO4'}}
if ligand_residues:
    resname_sel = ' or '.join([f'resname {{res}}' for res in ligand_residues])
    u.select_atoms(resname_sel).write('{output_file}')
else:
    # No ligands found, create empty file
    with open('{output_file}', 'w') as f:
        f.write('\\n')
"'''
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            raise Exception(f"MDAnalysis error: {result.stderr}")
        
        # If specific ligand residue name was provided, convert ATOM to HETATM
        if ligand_residue_name:
            convert_atom_to_hetatm(output_file)
        
        return True
    except Exception as e:
        print(f"Error in extract_ligands: {e}")
        return False

def convert_atom_to_hetatm(pdb_file):
    """Convert ATOM records to HETATM in PDB file"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        # Convert ATOM to HETATM
        converted_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                # Replace ATOM with HETATM
                converted_line = 'HETATM' + line[6:]
                converted_lines.append(converted_line)
            else:
                converted_lines.append(line)
        
        # Write back to file
        with open(pdb_file, 'w') as f:
            f.writelines(converted_lines)
        
        print(f"Converted ATOM records to HETATM in {pdb_file}")
        return True
    except Exception as e:
        print(f"Error converting ATOM to HETATM: {e}")
        return False

def correct_ligand_with_pymol(ligand_file, corrected_file):
    """Correct ligand using PyMOL"""
    ligand_path = os.path.abspath(ligand_file)
    corrected_path = os.path.abspath(corrected_file)
    if not os.path.isfile(ligand_path) or os.path.getsize(ligand_path) == 0:
        print("Ligand file missing or empty:", ligand_path)
        return False

    # Use PyMOL to add hydrogens and save corrected ligand
    cmd = f'pymol -cq {ligand_path} -d "h_add; save {corrected_path}; quit"'
    return run_command(cmd, f"Correcting ligand with PyMOL")

def remove_connect_records(pdb_file):
    """Remove CONNECT records from PDB file"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        # Filter out CONNECT records
        filtered_lines = [line for line in lines if not line.startswith('CONECT')]
        
        with open(pdb_file, 'w') as f:
            f.writelines(filtered_lines)
        
        print(f"Removed CONNECT records from {pdb_file}")
        return True
    except Exception as e:
        print(f"Error removing CONNECT records: {e}")
        return False

def merge_protein_and_ligand(protein_file, ligand_file, output_file):
    """Merge capped protein and corrected ligand with proper PDB formatting"""
    try:
        # Read protein file
        with open(protein_file, 'r') as f:
            protein_lines = f.readlines()
        
        # Read ligand file
        with open(ligand_file, 'r') as f:
            ligand_lines = f.readlines()
        
        # Process protein file: remove 'END' and add properly formatted 'TER'
        protein_processed = []
        last_atom_line = None
        for line in protein_lines:
            if line.strip() == 'END':
                # Create properly formatted TER card using the last atom's info
                if last_atom_line and last_atom_line.startswith('ATOM'):
                    # Extract atom number and residue info from last atom
                    atom_num = last_atom_line[6:11].strip()
                    res_name = last_atom_line[17:20].strip()
                    chain_id = last_atom_line[21:22].strip()
                    res_num = last_atom_line[22:26].strip()
                    ter_line = f"TER    {atom_num:>5}      {res_name} {chain_id}{res_num}\n"
                    protein_processed.append(ter_line)
                else:
                    protein_processed.append('TER\n')
            else:
                protein_processed.append(line)
                if line.startswith('ATOM'):
                    last_atom_line = line
        
        # Process ligand file: remove header info (CRYST, REMARK, etc.) and keep only ATOM/HETATM
        ligand_processed = []
        for line in ligand_lines:
            if line.startswith(('ATOM', 'HETATM')):
                ligand_processed.append(line)
        
        # Combine: protein + TER + ligand + END (no extra newline between TER and ligand)
        merged_content = ''.join(protein_processed) + ''.join(ligand_processed) + 'END\n'
        
        with open(output_file, 'w') as f:
            f.write(merged_content)
        
        return True
    except Exception as e:
        print(f"Error merging files: {str(e)}")
        return False

def prepare_structure(pdb_content, options, output_dir="output"):
    """Main function to prepare structure for AMBER simulation"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Define all file paths in output directory
        # Use completed structure if available, otherwise use original input
        complete_structure_file = os.path.join(output_dir, "0_complete_structure.pdb")
        original_input_file = os.path.join(output_dir, "0_original_input.pdb")
        
        if os.path.exists(complete_structure_file):
            input_file = complete_structure_file
            logger.info("Using completed structure (0_complete_structure.pdb) as input")
        else:
            input_file = original_input_file
            logger.info("Using original input (0_original_input.pdb) as input")
        
        user_chain_file = os.path.join(output_dir, "0_user_chain_selected.pdb")
        protein_file = os.path.join(output_dir, "1_protein_no_hydrogens.pdb")
        protein_capped_file = os.path.join(output_dir, "2_protein_with_caps.pdb")
        ligand_file = os.path.join(output_dir, "3_ligands_extracted.pdb")
        ligand_corrected_file = os.path.join(output_dir, "4_ligands_corrected.pdb")
        tleap_ready_file = os.path.join(output_dir, "tleap_ready.pdb")
        
        # Step 0: Save original input for reference (only if using original input)
        # If using completed structure, we don't overwrite it
        if input_file == original_input_file:
            print("Step 0: Saving original input...")
            with open(input_file, 'w') as f:
                f.write(pdb_content)
        else:
            # If using completed structure, read it instead of using pdb_content
            print("Step 0: Using completed structure as input...")
            with open(input_file, 'r') as f:
                pdb_content = f.read()
            # Also save a reference to original input if it doesn't exist
            if not os.path.exists(original_input_file):
                print("Step 0: Saving reference to original input...")
                with open(original_input_file, 'w') as f:
                    f.write(pdb_content)
        
        # Step 0.5: Extract user-selected chains and ligands
        selected_chains = options.get('selected_chains', [])
        selected_ligands = options.get('selected_ligands', [])
        
        if selected_chains:
            print(f"Step 0.5a: Extracting selected chains: {', '.join(selected_chains)}")
            if not extract_selected_chains(pdb_content, user_chain_file, selected_chains):
                raise Exception("Failed to extract selected chains")
        else:
            print("Step 0.5a: No chains selected, using original structure")
            shutil.copy2(input_file, user_chain_file)
        
        if selected_ligands:
            ligand_names = [f"{l.get('resn', '')}-{l.get('chain', '')}" for l in selected_ligands]
            print(f"Step 0.5b: Extracting selected ligands: {ligand_names}")
            if not extract_selected_ligands(pdb_content, ligand_file, selected_ligands):
                raise Exception("Failed to extract selected ligands")
        else:
            print("Step 0.5b: No ligands selected, creating empty ligand file")
            with open(ligand_file, 'w') as f:
                f.write('\n')
        
        # Step 1: Extract protein only (remove hydrogens) from user-selected chains
        print("Step 1: Extracting protein without hydrogens from selected chains...")
        # Read the user-selected chain file
        with open(user_chain_file, 'r') as f:
            chain_content = f.read()
        
        if not extract_protein_only(chain_content, protein_file):
            raise Exception("Failed to extract protein")
        
        # Step 2: Add capping groups (only if add_ace or add_nme is True)
        add_ace = options.get('add_ace', True)
        add_nme = options.get('add_nme', True)
        
        if add_ace or add_nme:
            print("Step 2: Adding ACE and NME capping groups...")
            if not add_capping_groups(protein_file, protein_capped_file):
                raise Exception("Failed to add capping groups")
        else:
            print("Step 2: Skipping capping groups (add_ace=False, add_nme=False)")
            print("Using protein without capping - copying to capped file")
            # Copy protein file to capped file (no capping)
            shutil.copy2(protein_file, protein_capped_file)
        
        # Step 3: Handle ligands (use pre-extracted ligand file)
        preserve_ligands = options.get('preserve_ligands', True)
        ligand_present = False
        
        if preserve_ligands:
            print("Step 3: Processing pre-extracted ligands...")
            
            # Check if ligand file has content (not just empty or newline)
            with open(ligand_file, 'r') as f:
                ligand_content = f.read().strip()
            
            if ligand_content and len(ligand_content) > 1:
                ligand_present = True
                print("Found pre-extracted ligands")
                
                # Correct ligand with PyMOL
                if not correct_ligand_with_pymol(ligand_file, ligand_corrected_file):
                    print("Error: Failed to process ligand")
                    return {
                        'error': 'Failed to process ligand with PyMOL',
                        'prepared_structure': '',
                        'original_atoms': 0,
                        'prepared_atoms': 0,
                        'removed_components': {},
                        'added_capping': {},
                        'preserved_ligands': 0,
                        'ligand_present': False
                    }
                
                # Merge protein and ligand
                if not merge_protein_and_ligand(protein_capped_file, ligand_corrected_file, tleap_ready_file):
                    raise Exception("Failed to merge protein and ligand")
            else:
                print("No ligands found in pre-extracted file, using protein only")
                # Copy protein file to tleap_ready
                shutil.copy2(protein_capped_file, tleap_ready_file)
        else:
            print("Step 3: Skipping ligand processing (preserve_ligands=False)")
            print("Using protein only - copying capped protein to tleap_ready")
            # Copy protein file to tleap_ready (protein only, no ligands)
            shutil.copy2(protein_capped_file, tleap_ready_file)
        
        # Remove CONNECT records from tleap_ready.pdb (PyMOL adds them)
        print("Removing CONNECT records from tleap_ready.pdb...")
        remove_connect_records(tleap_ready_file)
        
        # Read the final prepared structure
        with open(tleap_ready_file, 'r') as f:
            prepared_content = f.read()
            
            # Calculate statistics
            original_atoms = len([line for line in pdb_content.split('\n') if line.startswith('ATOM')])
            prepared_atoms = len([line for line in prepared_content.split('\n') if line.startswith('ATOM')])
            
            # Calculate removed components
            water_count = len([line for line in pdb_content.split('\n') if line.startswith('HETATM') and line[17:20].strip() in ['HOH', 'WAT', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'SPCE']])
            ion_count = len([line for line in pdb_content.split('\n') if line.startswith('HETATM') and line[17:20].strip() in ['NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 'CD', 'HG', 'PB', 'SR', 'BA', 'RB', 'CS', 'LI', 'F', 'BR', 'I', 'PO4', 'PO3', 'H2PO4', 'HPO4', 'H3PO4']])
            hydrogen_count = len([line for line in pdb_content.split('\n') if line.startswith('ATOM') and line[76:78].strip() == 'H'])
            
            # If not preserving ligands, count them as removed
            ligand_count = 0
            if not preserve_ligands and ligand_present:
                # Count ligands from the pre-extracted file
                with open(ligand_file, 'r') as f:
                    ligand_lines = [line for line in f if line.startswith('HETATM')]
                ligand_count = len(set(line[17:20].strip() for line in ligand_lines))
            
            removed_components = {
                'water': water_count,
                'ions': ion_count,
                'hydrogens': hydrogen_count,
                'ligands': ligand_count
            }
            
            # Calculate added capping groups (only if capping was performed)
            if add_ace or add_nme:
                # Count unique ACE and NME residues, not individual atoms
                ace_residues = set()
                nme_residues = set()
                
                for line in prepared_content.split('\n'):
                    if line.startswith('ATOM') and 'ACE' in line:
                        # Extract residue number to count unique ACE groups
                        res_num = line[22:26].strip()
                        ace_residues.add(res_num)
                    elif line.startswith('ATOM') and 'NME' in line:
                        # Extract residue number to count unique NME groups
                        res_num = line[22:26].strip()
                        nme_residues.add(res_num)
                
                added_capping = {
                    'ace_groups': len(ace_residues),
                    'nme_groups': len(nme_residues)
                }
            else:
                added_capping = {
                    'ace_groups': 0,
                    'nme_groups': 0
                }
            
            # Count preserved ligands from the pre-extracted file
            preserved_ligands = 0
            if ligand_present and preserve_ligands:
                with open(ligand_file, 'r') as f:
                    ligand_lines = [line for line in f if line.startswith('HETATM')]
                preserved_ligands = len(set(line[17:20].strip() for line in ligand_lines))
            
            result = {
                'prepared_structure': prepared_content,
                'original_atoms': original_atoms,
                'prepared_atoms': prepared_atoms,
                'removed_components': removed_components,
                'added_capping': added_capping,
                'preserved_ligands': preserved_ligands,
                'ligand_present': ligand_present,
                'separate_ligands': options.get('separate_ligands', False)
            }
            
            # If separate ligands is enabled and ligands are present, include ligand content
            if ligand_present and options.get('separate_ligands', False):
                with open(ligand_corrected_file, 'r') as f:
                    result['ligand_content'] = f.read()
            
            return result
        
    except Exception as e:
        return {
            'error': str(e),
            'prepared_structure': '',
            'original_atoms': 0,
            'prepared_atoms': 0,
            'removed_components': {},
            'added_capping': {},
            'preserved_ligands': 0,
            'ligand_present': False
        }

def parse_structure_info(pdb_content):
    """Parse structure information for display"""
    lines = pdb_content.split('\n')
    atom_count = 0
    chains = set()
    residues = set()
    water_molecules = 0
    ions = 0
    ligands = set()
    hetatoms = 0
    
    # Common water molecule names
    water_names = {'HOH', 'WAT', 'TIP3', 'TIP4', 'SPC', 'SPCE'}
    
    # Common ion names
    ion_names = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 'CD', 'HG', 'PB', 'SR', 'BA', 'RB', 'CS', 'LI', 'F', 'BR', 'I', 'PO4', 'PO3', 'H2PO4', 'HPO4', 'H3PO4','SO4'}
    
    # Common ligand indicators
    ligand_indicators = {'ATP', 'ADP', 'AMP', 'GDP', 'GTP', 'NAD', 'FAD', 'HEM', 'HEME', 'COA', 'SAM', 'PLP', 'THF', 'FMN', 'FAD', 'NADP', 'UDP', 'CDP', 'TDP', 'GDP', 'ADP', 'ATP'}

    for line in lines:
        if line.startswith('ATOM'):
            atom_count += 1
            chain_id = line[21:22].strip()
            if chain_id:
                chains.add(chain_id)
            
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            residues.add(f"{res_name}{res_num}")
        elif line.startswith('HETATM'):
            hetatoms += 1
            res_name = line[17:20].strip()
            
            if res_name in water_names:
                water_molecules += 1
            elif res_name in ion_names:
                ions += 1
            elif res_name in ligand_indicators:
                ligands.add(res_name)

    # Count unique water molecules
    unique_water_residues = set()
    for line in lines:
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            if res_name in water_names:
                unique_water_residues.add(f"{res_name}{res_num}")

    return {
        'atom_count': atom_count,
        'chains': list(chains),
        'residue_count': len(residues),
        'water_molecules': len(unique_water_residues),
        'ions': ions,
        'ligands': list(ligands),
        'hetatoms': hetatoms
    }

def test_structure_preparation():
    """Test function to verify structure preparation works correctly"""
    # Create a simple test PDB content
    test_pdb = """HEADER    TEST PROTEIN
ATOM      1  N   MET A   1      16.347  37.019  21.335  1.00 50.73           N  
ATOM      2  CA  MET A   1      15.737  37.120  20.027  1.00 45.30           C  
ATOM      3  C   MET A   1      15.955  35.698  19.546  1.00 41.78           C  
ATOM      4  O   MET A   1      16.847  35.123  20.123  1.00 40.15           O  
ATOM      5  CB  MET A   1      14.234  37.456  19.789  1.00 44.12           C  
ATOM      6  CG  MET A   1      13.456  36.123  19.234  1.00 43.45           C  
ATOM      7  SD  MET A   1      12.123  35.456  18.123  1.00 42.78           S  
ATOM      8  CE  MET A   1      11.456  34.123  17.456  1.00 42.11           C  
ATOM      9  N   ALA A   2      15.123  35.456  18.789  1.00 40.44           N  
ATOM     10  CA  ALA A   2      14.456  34.123  18.123  1.00 39.77           C  
ATOM     11  C   ALA A   2      13.123  33.456  17.456  1.00 39.10           C  
ATOM     12  O   ALA A   2      12.456  32.123  16.789  1.00 38.43           O  
ATOM     13  CB  ALA A   2      13.789  33.123  17.123  1.00 38.76           C  
ATOM     14  N   ALA A   3      12.789  32.456  16.123  1.00 38.09           N  
ATOM     15  CA  ALA A   3      11.456  31.789  15.456  1.00 37.42           C  
ATOM     16  C   ALA A   3      10.123  30.456  14.789  1.00 36.75           C  
ATOM     17  O   ALA A   3       9.456  29.123  14.123  1.00 36.08           O  
ATOM     18  CB  ALA A   3       9.789  29.456  13.456  1.00 35.41           C  
ATOM     19  OXT ALA A   3       8.123  28.789  13.456  1.00 35.74           O  
HETATM   20  O   HOH A   4      20.000  20.000  20.000  1.00 20.00           O  
HETATM   21  H1  HOH A   4      20.500  20.500  20.500  1.00 20.00           H  
HETATM   22  H2  HOH A   4      19.500  19.500  19.500  1.00 20.00           H  
HETATM   23  NA  NA  A   5      25.000  25.000  25.000  1.00 25.00          NA  
HETATM   24  CL  CL  A   6      30.000  30.000  30.000  1.00 30.00          CL  
HETATM    1  PG  GTP A 180      29.710  30.132  -5.989  1.00 52.48      A    P  
HETATM    2  O1G GTP A 180      29.197  28.937  -5.265  1.00 43.51      A    O  
HETATM    3  O2G GTP A 180      30.881  29.816  -6.827  1.00 63.11      A    O  
HETATM    4  O3G GTP A 180      30.013  31.278  -5.117  1.00 29.97      A    O  
HETATM    5  O3B GTP A 180      28.517  30.631  -6.995  1.00 23.23      A    O  
HETATM    6  PB  GTP A 180      27.017  31.171  -6.766  1.00 29.58      A    P  
HETATM    7  O1B GTP A 180      26.072  30.050  -6.958  1.00 17.62      A    O  
HETATM    8  O2B GTP A 180      26.960  31.913  -5.483  1.00 38.76      A    O  
HETATM    9  O3A GTP A 180      26.807  32.212  -7.961  1.00 13.12      A    O  
HETATM   10  PA  GTP A 180      26.277  33.726  -8.045  1.00 25.06      A    P  
HETATM   11  O1A GTP A 180      25.089  33.867  -7.187  1.00 44.06      A    O  
HETATM   12  O2A GTP A 180      27.427  34.635  -7.843  1.00 23.47      A    O  
HETATM   13  O5' GTP A 180      25.804  33.834  -9.555  1.00 42.05      A    O  
HETATM   14  C5' GTP A 180      26.615  33.475 -10.679  1.00 19.97      A    C  
HETATM   15  C4' GTP A 180      26.219  34.288 -11.894  1.00 14.90      A    C  
HETATM   16  O4' GTP A 180      24.826  34.017 -12.143  1.00 19.00      A    O  
HETATM   17  C3' GTP A 180      26.372  35.802 -11.724  1.00  4.96      A    C  
HETATM   18  O3' GTP A 180      26.880  36.347 -12.936  1.00 44.49      A    O  
HETATM   19  C2' GTP A 180      24.932  36.243 -11.481  1.00 17.12      A    C  
HETATM   20  O2' GTP A 180      24.719  37.581 -11.901  1.00 32.45      A    O  
HETATM   21  C1' GTP A 180      24.069  35.240 -12.240  1.00 16.17      A    C  
HETATM   22  N9  GTP A 180      22.724  35.005 -11.630  1.00 28.10      A    N  
HETATM   23  C8  GTP A 180      22.443  34.655 -10.325  1.00 27.05      A    C  
HETATM   24  N7  GTP A 180      21.168  34.483 -10.079  1.00 33.25      A    N  
HETATM   25  C5  GTP A 180      20.554  34.737 -11.307  1.00 26.23      A    C  
HETATM   26  C6  GTP A 180      19.183  34.712 -11.659  1.00 29.31      A    C  
HETATM   27  O6  GTP A 180      18.205  34.448 -10.957  1.00 40.80      A    O  
HETATM   28  N1  GTP A 180      19.000  35.036 -13.013  1.00 26.85      A    N  
HETATM   29  C2  GTP A 180      20.022  35.339 -13.903  1.00 28.70      A    C  
HETATM   30  N2  GTP A 180      19.627  35.619 -15.147  1.00 44.24      A    N  
HETATM   31  N3  GTP A 180      21.301  35.367 -13.569  1.00 21.67      A    N  
HETATM   32  C4  GTP A 180      21.489  35.054 -12.257  1.00 41.91      A    C  
END
"""
    
    options = {
        'remove_water': True,
        'remove_ions': True,
        'remove_hydrogens': True,
        'add_ace': True,
        'add_nme': True,
        'preserve_ligands': True,
        'separate_ligands': False,
        'fix_missing_atoms': False,
        'standardize_residues': False
    }
    
    print("Testing structure preparation...")
    result = prepare_structure(test_pdb, options, "output")
    
    print("\n=== STATISTICS ===")
    print(f"Original atoms: {result['original_atoms']}")
    print(f"Prepared atoms: {result['prepared_atoms']}")
    print(f"Removed: {result['removed_components']}")
    print(f"Added: {result['added_capping']}")
    print(f"Ligands: {result['preserved_ligands']}")
    print(f"Ligand present: {result['ligand_present']}")
    
    print(f"\nTest completed! Check 'output' folder for results:")
    print("- 1_protein_no_hydrogens.pdb (protein without hydrogens)")
    print("- 2_protein_with_caps.pdb (protein with ACE/NME caps)")
    print("- 3_ligands_extracted.pdb (extracted ligands, if any)")
    print("- 4_ligands_corrected.pdb (corrected ligands, if any)")
    print("- tleap_ready.pdb (final structure ready for tleap)")

if __name__ == "__main__":
    test_structure_preparation()