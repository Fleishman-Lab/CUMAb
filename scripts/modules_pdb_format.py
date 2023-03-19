"""
Modules for formatting pdb to be in the CUMAb format
Light chain Fv as chain A
heavy chain Fv as chain B
antigen chain as chain C
"""

import os
import sys
from typing import Dict, List, Tuple


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO, SearchIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.SeqUtils import seq1

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

CUMAb_dir = [x for x in sys.path if "CUMAb" in x.split("/")][0]
while CUMAb_dir.split("/")[-1] != "CUMAb":
    CUMAb_dir = CUMAb_dir.split("/" + CUMAb_dir.split("/")[-1])[0]

def fasta_from_pdb(pdb: str) -> None:
    name = pdb.split("/")[-1].split(".")[0]
    seq = ""
    structure = PDBParser().get_structure(name, pdb)
    for chain in structure.get_chains():
        count = 0
        het_prev = False
        for residue in chain:
            tags = residue.id
            if tags[0] != " ":
                if count != 0:
                    het_prev = True
            else:
                if het_prev == True:
                    sys.exit("HETATOM followed by ATOM within a chain")
                else:
                    seq += seq1(residue.get_resname())
    bio_seq_rec = SeqRecord(Seq.Seq(seq, generic_protein), id=name)
    SeqIO.write(bio_seq_rec, f'pdb_adjusting/{name}.fasta', "fasta")

def hmmscan(pdb: str) -> None:
    os.mkdir("pdb_adjusting/hmm")
    pdb_name = pdb.split("/")[-1].split(".")[0]
    hmmscan_cmd = f'hmmscan --domtblout '
    hmmscan_cmd +=  f'pdb_adjusting/hmm/{pdb_name}.dout'
    hmmscan_cmd += f' -o pdb_adjusting/hmm/{pdb_name}.out'
    hmmscan_cmd += f' {CUMAb_dir}/hmm_database/igg-chains_cdhit.hmm pdb_adjusting/{pdb_name}.fasta'
    os.system(hmmscan_cmd)

def extract_variable_region(pdb: str) -> None:
    pdb_name = pdb.split("/")[-1].split(".")[0]
    dout = f"pdb_adjusting/hmm/{pdb_name}.dout"
    dtbl = SearchIO.read(dout,'hmmscan3-domtab')
    light_max_score = 0
    heavy_max_score = 0
    kappa = ""
    for hmm in dtbl:
        dom = hmm.id.split("_")[0]
        for hit in hmm:
            score = hit.bitscore
            h = (dom == 'heavy') and (score > 80)
            l = dom in ['kappa' , 'lambda'] and (score > 80)
            if h and score > heavy_max_score:
                heavy_max_score = score
                heavy_coords = (hit.query_start,
                                        hit.query_end)
            if l and score > light_max_score:
                kappa = dom
                light_max_score = score
                light_coords = (hit.query_start,
                                        hit.query_end)
    
    heavy_range = range(heavy_coords[0], heavy_coords[1])
    light_range = range(light_coords[0], light_coords[1])
    assert light_coords[0] not in heavy_range
    assert light_coords[1] not in heavy_range

    with open(f"kappa_or_lambda.txt", "w") as f:
        f.write(kappa)

    fasta = f'pdb_adjusting/{pdb_name}.fasta'
    for record in SeqIO.parse(fasta, "fasta"):
        Ab_seq = str(record.seq)

    heavy_start = heavy_coords[0]
    heavy_end = heavy_coords[1]

    light_start = light_coords[0]
    light_end = light_coords[1]

    with open(f"pdb_adjusting/variable_coords.txt", "w") as f:
        f.write(f"{light_coords[0]} {light_coords[1]}\n")
        f.write(f"{heavy_coords[0]} {heavy_coords[1]}\n")

def cut_pdb(pdb: str, antigen_chain: str) -> None:
    cut_pdb_command = f"pymol -cq {CUMAb_dir}/scripts/pdbcutter.py -- {pdb}"
    cut_pdb_command += f" pdb_adjusting/variable_coords.txt {antigen_chain}"
    os.system(cut_pdb_command)

def combine_pdbs(pdb: str, antigen_chain: str) -> None:
    pdb_name = pdb.split("/")[-1].split(".")[0]
    light_pdb = f"pdb_adjusting/{pdb_name}_light.pdb"
    heavy_pdb = f"pdb_adjusting/{pdb_name}_heavy.pdb"
    antigen_pdb = f"pdb_adjusting/{pdb_name}_antigen.pdb"
    
    def add_lines(pdb:str):
        with open(pdb, "r") as f:
            for line in f:
                if "ATOM" in line or "TER" in line:
                    new.write(line)

    new_filename = f"{pdb_name}_CUMAb_format.pdb"
    with open(new_filename, "w") as new:
        add_lines(light_pdb)
        add_lines(heavy_pdb)
        if antigen_chain != "none":
            add_lines(antigen_pdb)
        new.write("END\n")

def format_pdb_file(pdb: str, antigen_chain: str) -> None:
    os.mkdir("pdb_adjusting")
    fasta_from_pdb(pdb)
    hmmscan(pdb)
    extract_variable_region(pdb)
    cut_pdb(pdb, antigen_chain)
    combine_pdbs(pdb, antigen_chain)
