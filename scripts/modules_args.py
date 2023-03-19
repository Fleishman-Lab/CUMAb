import json
import argparse
from typing import Dict, List, Tuple

def parse_arguments():
    parser = argparse.ArgumentParser(description="Takes input pdb of murine antibody, creates humanized designs")
    parser.add_argument('-pdb_file', nargs="?", type=str,help='Name of file of structure to use', required=True)
    parser.add_argument('-mode', nargs='?', type=str ,default="CDR",help='Mode to run. Options are CDR, SDR. Default is CDR')
    parser.add_argument("-antigen_chain", nargs="?", type=str, default="none", help="Name of antigen chain. Must be present in provided structure")
    parser.add_argument("-screens", nargs="+", type=str, default=[], help="Space separated list of regex expressions you would like to use as screens. CUMAb removes germlines that has these in the framework.")
    parser.add_argument("-origin_species", nargs="?", type=str, default="mouse", help="What species does your antibody come from? Options are mouse, human, rabbit. Default is mouse. Human and mouse are treated the same by CUMAb. Note that rabbit may be buggy and has not been validated experimentally")
    parser.add_argument("-res_to_fix", default = "none", nargs="+", type=str, help="Space separated list of residues to keep as mouse. Residues need to be formatted as residue number and then L or H for light/heavy chain. Ex: 5L or 6H. Residue number corresponds to the number in the sequence if the pdb is numbered from 1 with no gaps.")
    args = parser.parse_args()
    pdb_file = args.pdb_file
    mode = args.mode
    antigen_chain = args.antigen_chain
    screens = args.screens
    origin_species=args.origin_species
    res_to_fix = args.res_to_fix
    assert origin_species in ["human", "mouse", "rabbit"], "origin species must be either 'human', 'mouse', or 'rabbt'"
    assert mode in ["CDR", "SDR"], "mode must be either 'CDR' or 'SDR'"
    if not screens:
        screens = ["NG", "N[^P][ST]"]
    if mode == "SDR":
        assert antigen_chain, "SDR mode requires an antigen chain to be given"
    return  pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix

def write_config_JSON(mode: str, antigen_chain: str, 
                        screens: List, origin_species: str,
                        res_to_fix: List, pdb_file: str) -> None:
    args_dict = {"mode": mode, "antigen_chain": antigen_chain,
    "screens": screens, "origin_species":origin_species, 
    "res_to_fix": res_to_fix, "pdb_file": pdb_file}
    with open('config.json', 'w') as json_file:
        json.dump(args_dict, json_file)

def read_config():
    with open("config.json", "r") as json_file:
        config = json.load(json_file)
    mode = config["mode"]
    antigen_chain = config["antigen_chain"]
    screens = config["screens"]
    origin_species = config["origin_species"]
    res_to_fix = config["res_to_fix"]
    pdb_file = config["pdb_file"]
    return pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix
