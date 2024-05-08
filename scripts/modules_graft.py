"""
Modules for generating CDR-grafted sequences
"""

import os
import sys
import json
import re
from typing import Dict, List, Tuple
import pandas as pd

from Bio import SeqIO

CUMAb_dir = [x for x in sys.path if "CUMAb" in x.split("/")][0]
while CUMAb_dir.split("/")[-1] != "CUMAb":
    CUMAb_dir = CUMAb_dir.split("/" + CUMAb_dir.split("/")[-1])[0]

def fasta_database_parser() -> List:
    databases = ["IGHV", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]

    heavy_V_database_dict = {}
    heavy_J_database_dict = {}
    light_KV_database_dict = {}
    light_KJ_database_dict = {}
    light_LV_database_dict = {}
    light_LJ_database_dict = {}
    for database in databases:
        for record in SeqIO.parse(f"{CUMAb_dir}/IMGT_databases/{database}.fasta", "fasta"):
            name = record.name.split("|")[1]
            seq = str(record.seq)
            desc = record.description
            if "partial" not in desc and "rev" not in desc and desc.split("|")[3] == "F":
                if "H" in database:
                    if "V" in database:
                        if name.split("*")[0] not in heavy_V_database_dict:
                            if len(re.findall("C", seq)) <= 2:
                                heavy_V_database_dict[name.split("*")[0]] = seq
                    elif "J" in database:
                        if name.split("*")[0] not in heavy_J_database_dict:
                            heavy_J_database_dict[name.split("*")[0]] = seq
                elif "L" in database:
                    if "V" in database:
                        if name.split("*")[0] not in light_LV_database_dict:
                            if len(re.findall("C", seq)) <= 2:
                                light_LV_database_dict[name.split("*")[0]] = seq
                    elif "J" in database:
                        if name.split("*")[0] not in light_LJ_database_dict:
                            light_LJ_database_dict[name.split("*")[0]] = seq
                elif "K" in database:
                    if "V" in database:
                        if name.split("*")[0] not in light_KV_database_dict:
                            if len(re.findall("C", seq)) <= 2:
                                light_KV_database_dict[name.split("*")[0]] = seq
                    elif "J" in database:
                        if name.split("*")[0] not in light_KJ_database_dict:
                            light_KJ_database_dict[name.split("*")[0]] = seq
    return [heavy_V_database_dict, heavy_J_database_dict, light_KV_database_dict, light_KJ_database_dict, light_LV_database_dict, light_LJ_database_dict]

def make_full_chain(V_dict: Dict, J_dict: Dict) -> Dict:
    return_dict = {}
    for V in V_dict:
        for J in J_dict:
            key = V + "-" + J
            seq = V_dict[V] + J_dict[J]
            if len(re.findall("C", seq)) == 2:
                if key not in return_dict:
                    return_dict[key] = seq
    return return_dict

def parse_IMGT_databases() -> List:
    database_dicts = fasta_database_parser()
    heavy_V_database_dict = database_dicts[0]
    heavy_J_database_dict = database_dicts[1]
    light_KV_database_dict = database_dicts[2]
    light_KJ_database_dict = database_dicts[3]
    light_LV_database_dict = database_dicts[4]
    light_LJ_database_dict = database_dicts[5]
    heavy_dict = make_full_chain(heavy_V_database_dict, heavy_J_database_dict)
    kappa_dict = make_full_chain(light_KV_database_dict, light_KJ_database_dict)
    lambda_dict = make_full_chain(light_LV_database_dict, light_LJ_database_dict)
    return [heavy_dict, kappa_dict, lambda_dict]

def read_fasta(fasta:str) -> str:
    record = SeqIO.read(fasta, 'fasta')
    return str(record.seq)

def find_CDRs(type:str, sequence:str, species = "mouse") -> List:
    print("SEQUENCE:", sequence)
    CDRS = []
    assert type in ["light", "heavy"]
    if type == "light":
        L1 = re.search("C(.{8}.+?)W", sequence)[1]
        L2 = re.search(".{11}(.{10}.*?).{32}C", sequence.split(L1)[1])[1]
        L3 = re.search(".{24}.+?C(.+?)FG.G", sequence.split(L2)[1])[1]
        CDRS = [L1, L2, L3]
    if type == "heavy":
        H1 = re.search("C.{1}(.{11}.+?)W[^GNS]", sequence)[1]
        if species == "rabbit":
            if re.search(".{11}(.+?).{7}[KR][LIVFTA][TSIA]",
                           sequence.split(H1)[1]):
                H2 = re.search(".{11}(.+?).{7}[KR][LIVFTA][TSIA]",
                                sequence.split(H1)[1])
                if re.search("C(.+?)WG.G", sequence.split(H2)[1]):
                    H3 = re.search("C(.+?)WG.G", sequence.split(H2)[1])[1]
                elif re.search("C(.+?)WD.G", sequence.split(H2)[1]):
                    H3 = re.search("C(.+?)WD.G", sequence.split(H2)[1])[1]
                else:
                    sys.exit("Could not find CDR H3")
            else:
                sys.exit("Could not find CDR L1")
        else:
            H2 = re.search(".{11}(.+?).{36}C", sequence.split(H1)[1])[1]
            if re.search(".{36}C(.+?)WG.G", sequence.split(H2)[1]):
                H3 = re.search(".{36}C(.+?)WG.G", sequence.split(H2)[1])[1]
            elif re.search(".{36}C(.+?)WD.G", sequence.split(H2)[1]):
                H3 = re.search(".{36}C(.+?)WD.G", sequence.split(H2)[1])[1]
            else:
                sys.exit("Could not find CDR H3")
        CDRS = [H1, H2, H3]
    return CDRS

def check_CDR_lengths(mouse_light:List, mouse_heavy:List,
                        human_light:List, human_heavy:List) -> bool:
    to_return = True
    for i in range(3):
        if len(mouse_light[i]) != len(human_light[i]):
            to_return = False
    for j in range(2):
        if len(mouse_heavy[j]) != len(human_heavy[j]):
            to_return = False
    if len(human_heavy[2]) > len(mouse_heavy[2]):
        to_return = False
    return to_return

def graft_CDRs(sequence:str, mouse_CDRs:List,
                human_CDRs:List, mode="CDR") -> str:
    assert len(mouse_CDRs) == len(human_CDRs)
    if mode == "CDR":
        for i in range(len(mouse_CDRs)):
            sequence = sequence.replace(human_CDRs[i], mouse_CDRs[i])
    elif mode == "SDR":
        mouse_H3 = mouse_CDRs[2]
        human_H3 = human_CDRs[2]
        if len(mouse_H3) != len(human_H3):
            mouse_H3_list = list(mouse_H3)
            human_H3_list = list(human_H3)
            new_human_H3 = ""
            if len(mouse_H3) > len(human_H3):
                for res in range(2):
                    new_human_H3 += human_H3[res]
                for ins in range(2, 2+len(mouse_H3)-len(human_H3)):
                    new_human_H3 += mouse_H3[ins]
                for end_res in range(2, len(human_H3)):
                    new_human_H3 += human_H3[end_res]
            elif len(human_H3) > len(mouse_H3):
                if len(human_H3) - len(mouse_H3) == 1:
                    for res in range(0,7):
                        new_human_H3 += human_H3[res]
                    for res in range(8, len(human_H3)):
                        new_human_H3 += human_H3[res]
                elif len(human_H3) - len(mouse_H3) > 1:
                    return False
            sequence = sequence.replace(human_H3, new_human_H3)
    return sequence

def adjust_seq(type:str, mouse_seq:str, humanized_seq:str) -> str:
    assert type in ["light", "heavy"]
    mouse_before_C = len(re.search("(.+?)C", mouse_seq)[1])
    humanized_before_C = len(re.search("(.+?)C", humanized_seq)[1])
    before_C_dif = humanized_before_C - mouse_before_C
    if type == "light":
        mouse_end_len = len(re.search(".{90}FG.G(.+)", mouse_seq)[1])
        humanized_end_len = len(re.search(".{90}FG.G(.+)", humanized_seq)[1])
    if type == "heavy":
        if re.search(".{90}WG[^S]G(.+)", mouse_seq):
            mouse_end_len = len(re.search(".{90}WG[^S]G(.+)", mouse_seq)[1])
        elif re.search(".{90}WD[^S]G(.+)", mouse_seq):
            mouse_end_len = len(re.search(".{90}WD[^S]G(.+)", mouse_seq)[1])
        else:
            sys.exit("Cannot find mouse end length")
        humanized_end_len = len(re.search(".{90}WG[^S]G(.+)", humanized_seq)[1])
    end_dif = humanized_end_len - mouse_end_len
    sequence = humanized_seq[before_C_dif: len(humanized_seq) - end_dif]
    message = f"human seq is {sequence} and mouse seq is {mouse_seq}"
    assert len(sequence) == len(mouse_seq), message
    return sequence

def find_interface_locs(antigen_chain: str) -> str:
    interface_locs = ""
    if antigen_chain == "none":
        return interface_locs
    else:
        with open(f"interface_residues.txt") as f:
            for line in f:
                if "pymol" in line:
                    first_split = line.split("resi ")[1]
                    split = first_split.split("+")
                    for i in split:
                        if i != "+":
                            interface_locs += "{},".format(i.strip())
        interface_locs = interface_locs.strip(",")
        return interface_locs

def graft_interface_locs(human_seq:str, mouse_seq:str, locs:str) -> str:
    split = locs.split(",")
    list_human_seq = list(human_seq)
    if locs:
        for loc in split:
            int_loc = int(loc)
            list_human_seq[int_loc-1] = mouse_seq[int_loc-1]
    return "".join(list_human_seq)

def screen_motifs(light_seq:str, heavy_seq:str, screens:List) -> bool:
    value = False
    light_CDRs = find_CDRs("light", light_seq)
    heavy_CDRs = find_CDRs("heavy", heavy_seq)
    for screen in screens:
        r_string = r"{}".format(screen)
        light_count = 0
        heavy_count = 0
        for light_CDR in light_CDRs:
            light_count += len(re.findall(r_string, light_CDR))
        for heavy_CDR in heavy_CDRs:
            heavy_count += len(re.findall(r_string, heavy_CDR))
        if len(re.findall(r_string, light_seq)) > light_count:
            value = True
        if len(re.findall(r_string, heavy_seq)) > heavy_count:
            value = True
    return value

def graft_sequences(pdb_file: str, mode: str, antigen_chain: str, screens: List, origin_species: str, res_to_fix) -> None:
    name = pdb_file.split("/")[-1].split(".")[0]
    sequence_dicts = parse_IMGT_databases()
    heavy_dict = sequence_dicts[0]
    kappa_dict = sequence_dicts[1]
    lambda_dict = sequence_dicts[2]

    keys = []
    seqs = []
    with open(f"kappa_or_lambda.txt") as f:
        for line in f:
            if "kappa" in line:
                light_dict = kappa_dict
            elif "lambda" in line:
                light_dict = lambda_dict

    light_seq = read_fasta(f"pdb_adjusting/{name}_light.fasta")
    heavy_seq = read_fasta(f"pdb_adjusting/{name}_heavy.fasta")
    chain = light_seq + heavy_seq
    
    mouse_light_CDRs = find_CDRs("light", light_seq, origin_species)
    mouse_heavy_CDRs = find_CDRs("heavy", heavy_seq, origin_species)

    for heavy in heavy_dict:
        for light in light_dict:
            key = light + "-" + heavy
            germline_heavy_seq = heavy_dict[heavy]
            germline_light_seq = light_dict[light]
            human_light_CDRs = find_CDRs("light", germline_light_seq, origin_species)
            human_heavy_CDRs = find_CDRs("heavy", germline_heavy_seq, origin_species)
            if mode == "SDR":
                if check_CDR_lengths(mouse_light_CDRs, mouse_heavy_CDRs, human_light_CDRs, human_heavy_CDRs):
                    grafted_heavy_seq = graft_CDRs(germline_heavy_seq, mouse_heavy_CDRs, human_heavy_CDRs, "SDR")
                    grafted_light_seq = germline_light_seq
                    final_light_seq = adjust_seq("light", light_seq, grafted_light_seq)
                    final_heavy_seq = adjust_seq("heavy", heavy_seq, grafted_heavy_seq)
                    final_seq = final_light_seq + final_heavy_seq
                    interface_locs = find_interface_locs(antigen_chain)
                    if res_to_fix != "none":
                        for res in res_to_fix:
                            chain_id = res[-1]
                            num = int(res.split(chain_id)[0])
                            if chain_id == "H":
                                num += len(light_seq)
                            if interface_locs:
                                interface_locs += f",{num}"
                            else:
                                interface_locs += f"{num}"
                    final_seq = graft_interface_locs(final_seq, chain, interface_locs)
                    if not screen_motifs(final_seq[0:len(final_light_seq)], final_seq[len(final_light_seq):], screens):
                        if final_seq not in seqs:
                            seqs.append(final_seq)
                            keys.append(key)
            elif mode == "CDR":
                grafted_light_seq = graft_CDRs(germline_light_seq, mouse_light_CDRs, human_light_CDRs)
                grafted_heavy_seq = graft_CDRs(germline_heavy_seq, mouse_heavy_CDRs, human_heavy_CDRs)
                final_light_seq = adjust_seq("light", light_seq, grafted_light_seq)
                final_heavy_seq = adjust_seq("heavy", heavy_seq, grafted_heavy_seq)
                final_seq = final_light_seq + final_heavy_seq
                interface_locs = find_interface_locs(antigen_chain)
                if res_to_fix != "none":
                    for res in res_to_fix:
                        chain_id = res[-1]
                        num = int(res.split(chain_id)[0])
                        if chain_id == "H":
                            num += len(light_seq)
                        if interface_locs:
                            interface_locs += f",{num}"
                        else:
                            interface_locs += f"{num}"
                final_seq = graft_interface_locs(final_seq, chain, interface_locs)
                if not screen_motifs(final_seq[0:len(final_light_seq)], final_seq[len(final_light_seq):], screens):
                    if final_seq not in seqs:
                        seqs.append(final_seq)
                        keys.append(key)

    CDR_grafted_dict = {"Germline combination": keys, "Sequence": seqs}
    CDR_grafted_df = pd.DataFrame(CDR_grafted_dict)
    CDR_grafted_df.to_csv(f"{name}_grafted_sequences.csv", index=False)

def read_pdb(pdb_file: str) -> List:
    light_chain = ""
    heavy_chain = ""
    for record in SeqIO.parse(pdb_file, "pdb-atom"):
        print("record.annotations: ", record.annotations)
        print("record.seq: ", record.seq)
        if record.annotations["chain"] == "A":
            light_chain = str(record.seq)
        elif record.annotations["chain"] == "B":
            heavy_chain = str(record.seq)
    return [light_chain, heavy_chain]