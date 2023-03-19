import os
import shutil
import sys
from scripts.modules_args import parse_arguments, write_config_JSON
from scripts.modules_pdb_format import format_pdb_file
from scripts.modules_graft import find_interface_locs

CUMAb_dir = [x for x in sys.path if "CUMAb" in x.split("/")][0]
while CUMAb_dir.split("/")[-1] != "CUMAb":
    CUMAb_dir = CUMAb_dir.split("/" + CUMAb_dir.split("/")[-1])[0]

def main():
    #check required databases were downloaded and named properly
    for family in ["IGHV", "IGKV", "IGHJ", "IGKJ", "IGLV", "IGLJ"]:
        assert os.path.isfile(f"{CUMAb_dir}/IMGT_databases/{family}.fasta"), f"{family} database was not downloaded properly"
    
    #parse arguments
    pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix = parse_arguments()

    #find name of pdb_file, will be used throughout the run    
    pdb_name = pdb_file.split("/")[-1].split(".pdb")[0]
    
    #write configuration to a JSON file
    write_config_JSON(mode, antigen_chain, screens, origin_species, res_to_fix, pdb_file)

    #format pdb file: light chain variable as chain A and heavy chain variable as chain B, with antigen chain as C if present
    format_pdb_file(pdb_file, antigen_chain)
    name = pdb_file.split("/")[-1].split(".")[0]
    assert os.path.isfile(f"{name}_CUMAb_format.pdb"), "Pdb was not formatted properly"

if __name__ == '__main__':
    main()
