import requests
from collections import defaultdict
from textwrap import wrap
import re


def get_pdb_id_from_pdb_file(pdb_path):
    """
    Extract the 4-character PDB ID from a PDB file.

    By convention, PDB files have a line starting with 'HEADER' where
    columns 63–66 contain the PDB ID code.

    If that cannot be found, this function will raise a ValueError so
    that the pipeline fails loudly instead of silently doing the wrong thing.
    """
    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("HEADER") and len(line) >= 66:
                pdb_id = line[62:66].strip()
                if pdb_id:
                    return pdb_id.upper()

    raise ValueError(
        f"Could not determine PDB ID from file: {pdb_path}. "
        "Expected a 'HEADER' record with ID in columns 63–66."
    )

GRAPHQL_URL = "https://data.rcsb.org/graphql"

def detect_missing_residues(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()

    missing_by_chain = defaultdict(list)

    for line in response.text.splitlines():
        if line.startswith("REMARK 465"):
            parts = line.split()
            if len(parts) >= 5 and parts[2].isalpha():
                resname = parts[2]
                chain = parts[3]

                # Extract residue number (strip insertion code)
                match = re.match(r"(\d+)", parts[4])
                if match:
                    resnum = int(match.group(1))
                    missing_by_chain[chain].append((resname, resnum))

    return dict(missing_by_chain)

def get_chain_sequences(pdb_id):
    query = """
    query ChainSequences($pdb_id: String!) {
      entry(entry_id: $pdb_id) {
        polymer_entities {
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
          polymer_entity_instances {
            rcsb_polymer_entity_instance_container_identifiers {
              auth_asym_id
            }
          }
        }
      }
    }
    """

    r = requests.post(
        GRAPHQL_URL,
        json={"query": query, "variables": {"pdb_id": pdb_id}}
    )
    r.raise_for_status()

    chain_seqs = {}

    for entity in r.json()["data"]["entry"]["polymer_entities"]:
        seq = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
        for inst in entity["polymer_entity_instances"]:
            chain = inst[
                "rcsb_polymer_entity_instance_container_identifiers"
            ]["auth_asym_id"]
            chain_seqs[chain] = seq

    return chain_seqs

def write_fasta_for_missing_chains(pdb_id, chains_with_missing):
    filename = f"{pdb_id}_chains_with_missing.fasta"

    with open(filename, "w") as f:
        for chain, seq in chains_with_missing.items():
            f.write(f">{pdb_id.upper()}_{chain}\n")
            for line in wrap(seq, 60):
                f.write(line + "\n")

    print(f"Wrote FASTA: {filename}")

def run_esmfold(sequence):
    response = requests.post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/",
        data=sequence,
        timeout=300
    )
    response.raise_for_status()
    return response.text

def rebuild_pdb_with_esmfold(
    pdb_id,
    chains_to_replace,
    output_pdb=None,
    original_pdb_path=None,
):
    """
    pdb_id: str
        Original crystal structure object name (e.g. '3hhr')

    chains_to_replace: list[str]
        Chains that were missing residues and replaced by ESMFold
        Example: ['A', 'B', 'C']

    output_pdb: str, optional
        Output PDB filename.

    original_pdb_path: str, optional
        Path to the original PDB file that should be loaded into PyMOL
        as the reference object named `pdb_id`. If None, defaults to
        '../../output/0_original_input.pdb'.
    """

    from pymol import cmd

    # -----------------------------
    # 0. Load original PDB into PyMOL
    # -----------------------------
    if original_pdb_path is None:
        # Default to the pipeline output location
        original_pdb_path = "../output/0_original_input.pdb"

    print(f"Loading original PDB from {original_pdb_path} as object '{pdb_id}'")
    cmd.load(original_pdb_path, pdb_id)

    if output_pdb is None:
        output_pdb = f"{pdb_id}_rebuilt.pdb"

    # -----------------------------
    # 1. Align each ESMFold chain and fix chain IDs
    # -----------------------------
    for chain in chains_to_replace:
        esm_obj = f"{pdb_id}_chain_{chain}_esmfold"

        # Load the ESMFold-generated PDB for this chain as a PyMOL object
        esm_pdb_filename = f"{pdb_id}_chain_{chain}_esmfold.pdb"
        print(f"Loading ESMFold PDB {esm_pdb_filename} as object '{esm_obj}'")
        cmd.load(esm_pdb_filename, esm_obj)

        # ESMFold outputs everything as chain A by default.
        # Rename the chain in the loaded object to match the target chain ID.
        print(f"Renaming chain A -> {chain} in {esm_obj}")
        cmd.alter(esm_obj, f"chain='{chain}'")
        cmd.sort(esm_obj)  # Rebuild internal indices after alter

        align_cmd = (
            f"{esm_obj} and name CA",
            f"{pdb_id} and chain {chain} and name CA"
        )

        print(f"Aligning {esm_obj} to {pdb_id} chain {chain}")
        cmd.align(*align_cmd)

    # -----------------------------
    # 2. Build selection strings
    # -----------------------------
    chains_str = "+".join(chains_to_replace)

    esm_objs_str = " or ".join(
        f"{pdb_id}_chain_{chain}_esmfold"
        for chain in chains_to_replace
    )

    selection = (
        f"({pdb_id} and not chain {chains_str}) or "
        f"({esm_objs_str})"
    )

    # -----------------------------
    # 3. Create final model
    # -----------------------------
    cmd.select("final_model", selection)

    # -----------------------------
    # 4. Save rebuilt structure
    # -----------------------------
    cmd.save(output_pdb, "final_model")

    print(f"✅ Final rebuilt structure saved as: {output_pdb}")


if __name__ == "__main__":
    # Path to the original input PDB used by the pipeline
    original_pdb_path = "../output/0_original_input.pdb"

    # Automatically infer the PDB ID from the original PDB file,
    # instead of hard-coding it (e.g., '3hhr').
    pdb_id = get_pdb_id_from_pdb_file(original_pdb_path)
    print(f"Detected PDB ID from original file: {pdb_id}")

    # 1) Find missing residues for this structure
    missing = detect_missing_residues(pdb_id)
    chain_sequences = get_chain_sequences(pdb_id)

    chains_with_missing = {
        chain: chain_sequences[chain]
        for chain in missing
        if chain in chain_sequences
    }

    # 2) Write FASTA for chains with missing residues
    write_fasta_for_missing_chains(pdb_id, chains_with_missing)

    # 3) Run ESMFold for each chain and save results
    esmfold_results = {}
    chains_to_replace = []

    for chain, seq in chains_with_missing.items():
        print(f"Running ESMFold for chain {chain}")
        pdb_text = run_esmfold(seq)
        esmfold_results[chain] = pdb_text
        chains_to_replace.append(chain)
        # Save each chain
        with open(f"{pdb_id}_chain_{chain}_esmfold.pdb", "w") as f:
            f.write(pdb_text)

    # 4) Rebuild PDB in PyMOL using original structure and ESMFold chains
    rebuild_pdb_with_esmfold(
        pdb_id,
        chains_to_replace,
        original_pdb_path=original_pdb_path,
    )
    #rebuild_pdb_with_esmfold(
    #    pdb_id,
    #    ['B', 'C'],
    #    original_pdb_path=original_pdb_path,
    #                            )


