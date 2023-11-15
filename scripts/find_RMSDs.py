"""
takes as input a path to a CUMAb formatted pdb file,
and a path to directory of pdbs

outputs one csv with CDR RMSD for each pdb
"""

import glob
import pandas as pd
import sys
import os
from Bio import SeqIO
from modules_graft import read_pdb, find_CDRs
from modules_args import read_config
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

formatted_pdb = sys.argv[1]
path_to_pdbs = sys.argv[2].strip("/")

#load arguments from CUMAb_pdb_format.py
pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix = read_config()

#find name of pdb_file, will be used throughout the run    
pdb_name = pdb_file.split("/")[-1].split(".pdb")[0]

#read sequences from CUMAb formatted pdb
chain_seqs = read_pdb(formatted_pdb)
light = chain_seqs[0]
heavy = chain_seqs[1]

#find CDR positions
L_CDRs = find_CDRs("light", light)
H_CDRs = find_CDRs("heavy", heavy)

target_L1_first = light.find(L_CDRs[0]) + 1
target_L1_last = target_L1_first  + len(L_CDRs[0]) - 1

target_L2_first = light.find(L_CDRs[1]) + 1
target_L2_last = target_L2_first  + len(L_CDRs[1]) - 1

target_L3_first = light.find(L_CDRs[2]) + 1
target_L3_last = target_L3_first  + len(L_CDRs[2]) - 1

target_H1_first = heavy.find(H_CDRs[0]) + 1
target_H1_last = target_H1_first  + len(H_CDRs[0]) - 1

target_H2_first = heavy.find(H_CDRs[1]) + 1
target_H2_last = target_H2_first  + len(H_CDRs[1]) - 1

target_H3_first = heavy.find(H_CDRs[2]) + 1
target_H3_last = target_H3_first  + len(H_CDRs[2]) - 1

#define function to find RMSD of each CDR
def find_CDR_RMSDS(pdb_file:str):
    cmd.load(pdb_file, "target")
    cmd.cealign("target", "template")
    target_str =  "resi %s-%s "%(target_L1_first,target_L1_last)
    target_str += "and target and name C+O and chain A"
    template_str = "resi %s-%s "%(target_L1_first,target_L1_last)
    template_str += "and template and name C+O and chain A"
    cmd.select("target_L1_CO",target_str)
    cmd.select("model_L1_CO",template_str)
    L1_pair_fit_val = cmd.rms_cur("target_L1_CO","model_L1_CO")

    target_str =  "resi %s-%s "%(target_L2_first,target_L2_last)
    target_str += "and target and name C+O and chain A"
    template_str = "resi %s-%s "%(target_L2_first,target_L2_last)
    template_str += "and template and name C+O and chain A"
    cmd.select("target_L2_CO",target_str)
    cmd.select("model_L2_CO",template_str)
    L2_pair_fit_val = cmd.rms_cur("target_L2_CO","model_L2_CO")

    target_str =  "resi %s-%s "%(target_L3_first,target_L3_last)
    target_str += "and target and name C+O and chain A"
    template_str = "resi %s-%s "%(target_L3_first,target_L3_last)
    template_str += "and template and name C+O and chain A"
    cmd.select("target_L3_CO",target_str)
    cmd.select("model_L3_CO",template_str)
    L3_pair_fit_val = cmd.rms_cur("target_L3_CO","model_L3_CO")

    target_str =  "resi %s-%s "%(target_H1_first,target_H1_last)
    target_str += "and target and name C+O and chain B"
    template_str = "resi %s-%s "%(target_H1_first,target_H1_last)
    template_str += "and template and name C+O and chain B"

    cmd.select("target_H1_CO",target_str)
    cmd.select("model_H1_CO",template_str)
    H1_pair_fit_val = cmd.rms_cur("target_H1_CO","model_H1_CO")

    target_str =  "resi %s-%s "%(target_H2_first,target_H2_last)
    target_str += "and target and name C+O and chain B"
    template_str = "resi %s-%s "%(target_H2_first,target_H2_last)
    template_str += "and template and name C+O and chain B"
    cmd.select("target_H2_CO",target_str)
    cmd.select("model_H2_CO",template_str)
    H2_pair_fit_val = cmd.rms_cur("target_H2_CO","model_H2_CO")

    target_str =  "resi %s-%s "%(target_H3_first,target_H3_last)
    target_str += "and target and name C+O and chain B"
    template_str = "resi %s-%s "%(target_H3_first,target_H3_last)
    template_str += "and template and name C+O and chain B"
    cmd.select("target_H3_CO",target_str)
    cmd.select("model_H3_CO",template_str)
    H3_pair_fit_val = cmd.rms_cur("target_H3_CO","model_H3_CO")
    cmd.delete("target")
    to_return = [L1_pair_fit_val, L2_pair_fit_val, L3_pair_fit_val]
    to_return += [H1_pair_fit_val, H2_pair_fit_val, H3_pair_fit_val]
    return to_return

files = glob.glob(f"{path_to_pdbs}/*")
pdbs = [x for x in files if ".pdb" in x]
L1_rmsd = []
L2_rmsd = []
L3_rmsd = []
H1_rmsd = []
H2_rmsd = []
H3_rmsd = []
names = []
cmd.load(relax_pdb, "template")
for pdb in pdbs:
    name = pdb.split("/")[-1].split(".pdb")[0]
    names.append(name)
    if os.path.getsize(file) > 0:
        try:
            rmsds = find_CDR_RMSDS(pdb_file)
            L1_rmsd.append(rmsds[0])
            L2_rmsd.append(rmsds[1])
            L3_rmsd.append(rmsds[2])
            H1_rmsd.append(rmsds[3])
            H2_rmsd.append(rmsds[4])
            H3_rmsd.append(rmsds[5])
        except:
            L1_rmsd.append(0)
            L2_rmsd.append(0)
            L3_rmsd.append(0)
            H1_rmsd.append(0)
            H2_rmsd.append(0)
            H3_rmsd.append(0)
final_dict = {"Sequence": names, "L1 RMSD": L1_rmsd, "L2 RMSD": L2_rmsd,
              "L3 RMSD": L3_rmsd, "H1 RMSD": H1_rmsd, "H2 RMSD": H2_rmsd, 
              "H3": H3_rmsd}
final_df = pd.DataFrame(final_dict)
final_df.to_csv(f"{pdb_name}_RMSDS.csv", index=False)
