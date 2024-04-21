import os
import shutil
import sys
from scripts.modules_graft import graft_sequences, find_interface_locs
from scripts.modules_args import read_config

def main():
	assert os.path.isfile("config.json"), "Missing config file. Please run CUMAb_pdb_format.py"

	#load arguments from CUMAb_pdb_format.py
	pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix = read_config()

	#create csv of grafted sequences
	graft_sequences(pdb_file, mode, antigen_chain, screens, origin_species, res_to_fix)


	name = pdb_file.split("/")[-1].split(".")[0]

	#write example command line for relax
	interface_locs = find_interface_locs(antigen_chain)
	relax_command_line = "{path_to_rosetta_exec} -database {path_to_rosetta_database} -parser:protocol {path_to_CUMAb_dir}/xmls/Relax.xml"
	relax_command_line += f" -s {name}_CUMAb_format.pdb"
	relax_command_line += " @{path_to_CUMAb_dir}/flags"
	relax_command_line += f" -parser:script_vars interface_residues={interface_locs}"
	with open("relax_command_line_example", "w") as f:
		f.write(relax_command_line)

	#write example command line for threading
	thread_command_line = "{path_to_rosetta_exec} -database {path_to_rosetta_database} -parser:protocol {path_to_CUMAb_dir}/xmls/CUMAb.xml -s {path_to_relaxed_pdb}"
	thread_command_line += " @{path_to_CUMAb_dir}/flags -parser:script_vars Sequence={humanized_seq} -overwrite"
	thread_command_line += f" -parser:script_vars interface_residues={interface_locs}"
	with open("thread_command_line_example", "w") as f:
		f.write(thread_command_line)

	#delete unused files
	assert os.path.isfile(f"{name}_grafted_sequences.csv"), "Error in CDR grafting"
	assert os.path.isfile("relax_command_line_example"), "Error in writing example relax command line"
	assert os.path.isfile("thread_command_line_example"), "Error in writing example thread command line"
	files_to_delete = [f"{name}_CUMAb_format_0001.pdb", "interface_residues.txt", "kappa_or_lambda.txt", "score.sc"]
	for file in files_to_delete:
		if os.path.isfile(file):
			os.remove(file)
	if os.path.isdir("pdb_adjusting"):
		shutil.rmtree("pdb_adjusting")

if __name__ == '__main__':
    main()