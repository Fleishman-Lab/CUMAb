# CUMAb
Licensed under the Non-Profit Open Software License version 3.0.

This repository contains the scripts and xmls needed to run the CUMAb tool for antibody humanization locally. The method is described in detail in our paper published in Nature Biomedical Engineering ([Link to paper](https://rdcu.be/diKmJ)). Please note that we also provide CUMAb as a webserver free to academics ([Link to webserver](https://CUMAb.weizmann.ac.il)) and that for most uses it is simpler to use the webserver. This repository is intended to allow users to generate the humanized sequences that CUMAb uses and provides an example command line to run the Rosetta calculations; however, it is up to the user to decide how to run the command line for each humanized sequence efficiently. Questions, comments, and suggestions can be sent to ariel.tennenhouse@weizmann.ac.il. 

## Citations
Please cite our manuscript as well as IMGT 
- Giudicelli, V., Chaume, D. & Lefranc, M.-P. IMGT/GENE-DB: a comprehensive database for human and mouse immunoglobulin and T cell receptor genes. Nucleic Acids Res. 33, D256–61 (2005).
- Tennenhouse, A., Khmelnitsky, L., Khalaila, R. et al. Computational optimization of antibody humanness and stability by systematic energy-based ranking. Nat. Biomed. Eng (2023). https://doi.org/10.1038/s41551-023-01079-1

## Installation
You will need to either have Rosetta installed or install it from http://www.rosettacommons.org. CUMAb uses git version d9d4d5dd3fd516db1ad41b302d147ca0ccd78abd
<br>
<br>You will need to download the IMGT databases of antibody germline sequences. Please save them under a folder called "IMGT_databases" in 6 separate files ("IGHV.fasta", "IGHJ.fasta", "IGKV.fasta", "IGKJ.fasta", "IGLJ.fasta", "IGKJ.fasta"). 
<br>They are available from the following urls:
  - [IGHV](https://www.imgt.org/genedb/GENElect?query=7.6+IGHV&species=Homo+sapiens)
  - [IGHJ](https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=Homo+sapiens)
  - [IGKV](https://www.imgt.org/genedb/GENElect?query=7.6+IGKV&species=Homo+sapiens)
  - [IGKJ](https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=Homo+sapiens)
  - [IGLV](https://www.imgt.org/genedb/GENElect?query=7.6+IGLV&species=Homo+sapiens)
  - [IGLJ](https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=Homo+sapiens)
<br>
You need to install the CUMAb conda environment. To do so, run:

```
conda env create -f CUMAb_environment.yml -n CUMAb
conda activate CUMAb
```

## Running CUMAb
- CUMAb takes as an input only a pdb file of the antibody you wish to humanize
  - If you do not have a structure of the antibody, we reccomend to use [ABodyBuilder2](https://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/abodybuilder2/) to predict the structure
- The final output includes:
  - a pdb file formatted as necessary for CUMAb
  - a csv file containing the humanized sequences
  - example command lines to run the initial relax on the pdb as well as the threading of the humanized sequences

### Step 1: Formatting the pdb file
- Create a new directory and move to that directory
- Place the pdb file of the antibody you wish to humanize in this directory
- Run the following command:
```
python {path_to_CUMAb_dir}/CUMAb_pdb_format.py -pdb_file {path_to_PDB_file}
```
- Arguments:
  - -pdb_file: Required. Path to pdb file of antibody you want to humanize. This should be the only thing in the new directory you made for the CUMAb run
  - -mode: Not required. Do you want to graft the sequences of the entire CDRs (CUMAb definition) or only SDRs? If you want to use SDR grafting, antigen chain must be provied in the pdb and as an argument here. Default is CDR and the other option is SDR
  - -antigen_chain: Not required. Does your pdb contain the antigen as well? Must be a single letter corresponding to the chain of the antigen in the pdb file. Default is none
  - -screens: Not required. Space separated list of regular expressions that sequences will be exlcuded if they contain outside of the CDRs. Default is NG N[^P][ST]
  - -origin_species: Not required. What species does your antibody originate from? Options are human, mouse, or rabbit. Human and mouse are treated the same whereas rabbit has some slight differences. Default is mouse
  - -res_to_fix: Not required. Space separated list of residues which should be kept the same identity as in the parental antibody. Must be formatted as residue number followed by either L or H for light or heavy chain. Residue number must be numbered from start of pdb with no gaps.

### Step 1b: Finding residues in antibody-antigen interface
- **Do not run if your pdb file does not contain an antigen**
- If your pdb file also contains the antigen of interest, run the following command to find residues in the antibody-antigen interface:
```
{path_to_rosetta_exec} -database {path_to_rosetta_database} -s {CUMAb_formatted_pdb} -parser:protocol {path_to_CUMAb_dir}/xmls/Interface.xml -overwrite | grep protocols.protein_interface_design.filters.DesignableResiduesFilter > interface_residues.txt 
```
### Step 2: Creating humanized sequences
From the same directory, fun the following command:
```
python {path_to_CUMAb_dir}/CUMAb_graft_sequences.py
```

### Step 3: Relax structure of parental antibody
We reccomend to relax the starting pdb structure before threading the humanized sequences onto it. An example command line for running the relax was created under the file name "relax_command_line_example". Please note that in our protocol for CUMAb we run the initial relax 15 times and take the lowest scoring one to use for the threading. 

### Step 4: Threading humanized sequences onto structure of parental antibody
An example command line for running the threading can be found under the file name "thread_command_example". You will need to run this command for each of the humanized sequences created in part 2. 
