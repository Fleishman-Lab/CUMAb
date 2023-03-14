# CUMAb
Scripts and xmls needed to run the CUMAb tool for antibody humanization as described in Nature Biomedical Engineering 2023. Please note that we also provide CUMAb as a webserver free to academics (www.CUMAb.weizmann.ac.il) and that for most uses it is simpler to use the webserver. This repository is intended to allow users to generate the humanized sequences and a command line for running each of them that CUMAb would use; however, it is up to the user to decide how to run these command lines efficiently. 

## Installation
You will need to either have Rosetta installed or install it from http://www.rosettacommons.org.
<br>
<br>You will need to download the IMGT databases of antibody germline sequences. Please save them under a folder called "IMGT_databases" in 6 separate files ("IGHV.fasta", "IGHJ.fasta", "IGKV.fasta", "IGKJ.fasta", "IGLJ.fasta", "IGKJ.fasta"). 
<br>They are available from the following urls:
  - IGHV: https://www.imgt.org/genedb/GENElect?query=7.6+IGHV&species=Homo+sapiens
  - IGHJ: https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=Homo+sapiens
  - IGKV: https://www.imgt.org/genedb/GENElect?query=7.6+IGKV&species=Homo+sapiens
  - IGKJ: https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=Homo+sapiens
  - IGLV: https://www.imgt.org/genedb/GENElect?query=7.6+IGLV&species=Homo+sapiens
  - IGLJ: https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=Homo+sapiens
 <br>
 <br> You need to install the CUMAb conda environment. To do so, run:
 '''
 conda create --name CUMAb --file CUMAb_spec.txt
 conda activate CUMAb
 '''
