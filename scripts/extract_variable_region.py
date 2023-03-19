"""
script that takes as input the .dout file from hmmscan and the fasta of the
full pdb as inputted by the user
the script searches the .dout file for the highest aligning chains
and finds the coordinates of the variable region
{Ab} = the name of .dout file
coordinates are outputted to {Ab}/variable_coords.txt
also finds if antibody is kappa or lambda and writes to {Ab}/korl.txt
Note: often cuts off variable region one residue on either side of
each chain
"""

from Bio import SeqIO, SearchIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import sys


dout = sys.argv[1]
fasta = sys.argv[2]

name = dout.split("/")[-1].split(".")[0]
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

with open(f".kappa_or_lambda.txt", "w") as f:
    f.write(kappa)

for record in SeqIO.parse(fasta, "fasta"):
    Ab_seq = str(record.seq)

heavy_start = heavy_coords[0]
heavy_end = heavy_coords[1]

light_start = light_coords[0]
light_end = light_coords[1]

with open(f"./variable_coords.txt", "w") as f:
    f.write(f"{light_coords[0]} {light_coords[1]}\n")
    f.write(f"{heavy_coords[0]} {heavy_coords[1]}\n")

