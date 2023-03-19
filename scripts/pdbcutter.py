"""
takes as input a pdb file, the file with coordinates of light and heavy chain
variable region from extract_variable region, and the antigen chain if there is
one

outputs a few files:
    fasta with light chain sequence
    fasta with heavy chain sequence
    pdb with light chain
    pdb with heavy chain
    pdb with antigen chain if there is
"""
import sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import re
import shutil

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

pdb = sys.argv[-3]
coords = sys.argv[-2]
antigen_chain = sys.argv[-1]

name = pdb.split("/")[-1].split(".")[0]
fasta = f"pdb_adjusting/{name}.fasta"
possible_chains = ["D", "E", "F", "G", "H", "I"]

def remove_hetero_atoms():
    cmd.select("not_hetatm", "not hetatm")
    cmd.save(pdb, "not_hetatm")

def rename_other_chains(chain:str, light:str, heavy:str, ant:str):
    if heavy != chain and light != chain and ant != chain:
        for poss in possible_chains:
            if heavy_chain != poss and light_chain != poss:
                if antigen_chain != poss:
                    cmd.alter(f"{name} and chain {chain}", f"chain='{poss}'")
                    break

def renumber_pdb_chain(chain:str):
    cmd.select("all_", f"name ca and chain {chain}")
    all = {"resnums":[]}
    cmd.iterate("all_",  "resnums.append(resi)", space=all)
    residues = list(all.get("resnums"))
    count = 1
    for i in residues:
        if not re.search("[a-zA-z]", str(i)):
            alter_str = "resi = '" + str(10000 + int(i)) + "'"
            cmd.alter(f"chain {chain} and resi {i}", alter_str)
    for i in residues:
        alter_str = "resi = '" + str(count) + "'"
        if not re.search("[a-zA-z]", str(i)):
            cmd.alter(f"chain {chain} and resi {int(i) + 10000}", alter_str)
        else:
            cmd.alter(f"chain {chain} and resi {i}", alter_str)
        count += 1

def check_interface(chain1:str, chain2:str):
    if not chain1 or not chain2:
        return True
    cmd.select("chain1", f"chain {chain1}")
    cmd.select("chain2", f"chain {chain2}")
    cmd.select("interface", f"chain1 and (chain2 around 10)")
    interface_list = {"resnums":[]}
    cmd.iterate("interface",  "resnums.append(resi)", space=interface_list)
    residues = interface_list.get("resnums")
    if residues:
        return True
    else:
        return False

def cut_str(the_string: str, start:int, end:int) -> str:
    return the_string[start:end]

count = 0
with open(coords, "r") as f:
    for line in f:
        split = line.split(" ")
        if count == 0:
            light_coords = [int(split[0]), int(split[1])]
        else:
            heavy_coords = [int(split[0]), int(split[1])]
        count += 1

korl = ""
with open("kappa_or_lambda.txt") as f:
    for line in f:
        if "kappa" in line:
            korl = "kappa"
        elif "lambda" in line:
            korl = "lambda"

for record in SeqIO.parse(fasta, "fasta"):
    Ab_seq = str(record.seq)
light_seq = Ab_seq[light_coords[0] : light_coords[1]+1]
light_len = len(light_seq)
heavy_seq = Ab_seq[heavy_coords[0] : heavy_coords[1]+1]
heavy_len = len(heavy_seq)
heavy_chain = ""
light_chain = ""
cmd.load(pdb,name)
remove_hetero_atoms()
cmd.delete(name)
cmd.load(pdb,name)
all_chains = []

for record in SeqIO.parse(pdb, "pdb-atom"):
    all_chains.append(record.annotations["chain"])
    renumber_pdb_chain(record.annotations["chain"])
    cmd.save(pdb,name)
for record in SeqIO.parse(pdb, "pdb-atom"):
    if heavy_seq in record.seq:
        if not heavy_chain:
            if not light_chain or check_interface(light_chain,
            record.annotations["chain"]):
                heavy_chain = record.annotations["chain"]
                heavy_start = str(record.seq).find(heavy_seq) + 1
                final_heavy_seq = cut_str(str(record.seq), heavy_start - 1,
                                    heavy_start + heavy_len)
                before_C_len = len(re.search("(.+?)C", final_heavy_seq)[1])
                while before_C_len < 21:
                    if heavy_start != 1:
                        heavy_start -= 1
                        heavy_len += 1
                        final_heavy_seq = cut_str(str(record.seq),
                                        heavy_start - 1,
                                         heavy_start + heavy_len)
                        before_C_len = len(re.search("(.+?)C",
                                            final_heavy_seq)[1])
                    else:
                        break
                if re.search("WG.G(.+$)", final_heavy_seq):
                    end_len = len(re.search("WG.G(.+$)",
                                            final_heavy_seq)[1])
                elif re.search("WD.G(.+$)", final_heavy_seq):
                    end_len = len(re.search("WD.G(.+$)",
                                            final_heavy_seq)[1])
                else:
                    sys.exit("Could not find end length in pdbcutter")

                while end_len < 7:
                    if len(record.seq) >= heavy_len + 1:
                        heavy_len += 1
                        final_heavy_seq = cut_str(str(record.seq),
                                                    heavy_start - 1,
                                                    heavy_start + heavy_len)
                        if re.search("WG.G(.+$)", final_heavy_seq):
                            end_len = len(re.search("WG.G(.+$)",
                                        final_heavy_seq)[1])
                        elif re.search("WD.G(.+$)", final_heavy_seq):
                            end_len = len(re.search("WD.G(.+$)",
                                                    final_heavy_seq)[1])
                        else:
                            sys.exit("Could not find end length in pdbcutter")
                    else:
                        break
                if end_len > 7:
                    heavy_len -= 1
                    final_heavy_seq = cut_str(str(record.seq),
                                            heavy_start - 1,
                                            heavy_start + heavy_len)
                heavy_start += record.annotations["start"] - 1
    elif heavy_seq[:-1] in record.seq:
        if not heavy_chain:
            if not light_chain or check_interface(light_chain,
                                            record.annotations["chain"]):
                heavy_chain = record.annotations["chain"]
                heavy_start = str(record.seq).find(heavy_seq[:-1]) + 1
                final_heavy_seq = cut_str(str(record.seq),heavy_start - 1,
                                            heavy_start + heavy_len - 1)
                heavy_start += record.annotations["start"] - 1
                before_C_len = len(re.search("(.+?)C", final_heavy_seq)[1])
                while before_C_len < 21:
                    if heavy_start != 1:
                        heavy_start -= 1
                        heavy_len += 1
                        final_heavy_seq = cut_str(str(record.seq),
                                            heavy_start - 1,
                                            heavy_start + heavy_len - 1)
                        before_C_len = len(re.search("(.+?)C",
                                            final_heavy_seq)[1])
                    else:
                        break
                if re.search("WG.G(.+$)", final_heavy_seq):
                    end_len = len(re.search("WG.G(.+$)",
                                            final_heavy_seq)[1])
                elif re.search("WD.G(.+$)", final_heavy_seq):
                    end_len = len(re.search("WD.G(.+$)",
                                            final_heavy_seq)[1])
                else:
                    sys.exit("Could not find end length in pdbcutter")
                while end_len > 7:
                    heavy_len -= 1
                    final_heavy_seq = cut_str(str(record.seq),
                                                heavy_start - 1,
                                                heavy_start + heavy_len)
                    if re.search("WG.G(.+$)", final_heavy_seq):
                        end_len = len(re.search("WG.G(.+$)",
                                    final_heavy_seq)[1])
                    elif re.search("WD.G(.+$)", final_heavy_seq):
                        end_len = len(re.search("WD.G(.+$)",
                                                final_heavy_seq)[1])
                    else:
                        sys.exit("Could not find end length in pdbcutter")
                heavy_start += record.annotations["start"] - 1

    if light_seq in record.seq:
        if not light_chain:
            if not heavy_chain or check_interface(heavy_chain,
                                            record.annotations["chain"]):
                light_chain = record.annotations["chain"]
                light_start = str(record.seq).find(light_seq) + 1
                if light_start != 1:
                    light_start -= 1
                    light_len += 1
                final_light_seq = cut_str(str(record.seq), light_start - 1,
                                            light_start + light_len)
                end_len = len(re.search("FG.G(.+$)", final_light_seq)[1])
                before_C_len = len(re.search("(.+?)C", final_light_seq)[1])
                if korl == "kappa":
                    while before_C_len < 22:
                        if light_start != 1:
                            light_len += 1
                            light_start -= 1
                            final_light_seq = cut_str(str(record.seq),
                                                    light_start - 1,
                                                    light_start + light_len)
                            before_C_len = len(re.search("(.+?)C",
                                                final_light_seq)[1])
                        else:
                            break

                elif korl == "lambda":
                    while before_C_len < 21:
                        if light_start != 1:
                            light_len += 1
                            light_start -= 1
                            final_light_seq = cut_str(str(record.seq),
                                                    light_start - 1,
                                                    light_start + light_len)
                            before_C_len = len(re.search("(.+?)C",
                                                final_light_seq)[1])
                        else:
                            break
                before_C_len = len(re.search("(.+?)C", final_light_seq)[1])
                while end_len < 6:
                    if len(record.seq) >= light_len + 1:
                        light_len += 1
                        final_light_seq = cut_str(str(record.seq),
                                                light_start - 1,
                                                light_start + light_len)
                        end_len = len(re.search("FG.G(.+$)",
                                        final_light_seq)[1])
                    else:
                        break
                while end_len > 6:
                    light_len -= 1
                    final_light_seq = cut_str(str(record.seq),
                                                light_start - 1,
                                                light_start + light_len)
                    end_len = len(re.search("FG.G(.+$)",
                                            final_light_seq)[1])
                light_start += record.annotations["start"] - 1

    elif light_seq[:-1] in record.seq:
        if not light_chain:
            light_chain = record.annotations["chain"]
            light_start = str(record.seq).find(light_seq[:-1]) + 1
            if light_start != 1:
                light_start -= 1
                light_len += 1
            final_light_seq = cut_str(str(record,seq), light_start - 1,
                                        light_start + light_len - 1)
            before_C_len = len(re.search("(.+?)C", final_light_seq)[1])
            count = 1
            if korl == "kappa":
                while before_C_len < 22:
                    if light_start != 1:
                        light_len += 1
                        light_start -= 1
                        final_light_seq = cut_str(str(record.seq),
                                                light_start - 1,
                                                light_start + light_len - 1)
                        before_C_len = len(re.search("(.+?)C",
                                            final_light_seq)[1])
                    else:
                        break
            elif korl == "lambda":
                while before_C_len < 21:
                    if light_start != 1:
                        light_len += 1
                        light_start -= 1
                        final_light_seq = cut_str(str(record.seq),
                                                light_start - 1,
                                                light_start + light_len - 1)
                        before_C_len = len(re.search("(.+?)C",
                                            final_light_seq)[1])
                    else:
                        break
            end_len = len(re.search("FG.G(.+$)", final_light_seq)[1])
            while end_len > 6:
                light_len -= 1
                final_light_seq = cut_str(str(record.seq),
                                            light_start - 1,
                                            light_start + light_len)
                end_len = len(re.search("FG.G(.+$)", final_light_seq)[1])
            light_start += record.annotations["start"] - 1

final_light_seq = final_light_seq.replace("X", "")
final_heavy_seq = final_heavy_seq.replace("X", "")
light_seq_bio = SeqRecord(Seq.Seq(final_light_seq, generic_protein), id=name)
heavy_seq_bio = SeqRecord(Seq.Seq(final_heavy_seq, generic_protein), id=name)
SeqIO.write(light_seq_bio, f'pdb_adjusting/{name}_light.fasta', "fasta")
SeqIO.write(heavy_seq_bio, f'pdb_adjusting/{name}_heavy.fasta', "fasta")
cmd.delete(name)
cmd.load(pdb,name)
light_len = len(final_light_seq)
heavy_len = len(final_heavy_seq)

rename_other_chains("A", light_chain, heavy_chain, antigen_chain)
rename_other_chains("B", light_chain, heavy_chain, antigen_chain)
if antigen_chain != "none":
    rename_other_chains("C", light_chain, heavy_chain, antigen_chain)
    changed = ""
    if antigen_chain == "A" or antigen_chain == "B":
        for ch in possible_chains:
            if ch != light_chain and ch != heavy_chain:
                cmd.alter(f"{name} and chain {antigen_chain}",
                            f"chain='{ch}'")
                changed = ch
                break
if light_chain != "B":
    cmd.alter(f"{name} and chain {heavy_chain} and resi {1}-{heavy_len}",
                "chain='B'")
    cmd.alter(f"{name} and chain {heavy_chain} and resi {1+heavy_len}-",
                "chain='Z'")
    cmd.alter(f"{name} and chain {light_chain} and resi {1}-{light_len}",
                "chain='A'")
    cmd.alter(f"{name} and chain {light_chain} and resi {1+light_len}-",
                "chain='Z'")
elif heavy_chain != "A":
    cmd.alter(f"{name} and chain {light_chain} and resi {1}-{light_len}",
                "chain='A'")
    cmd.alter(f"{name} and chain {light_chain} and resi {1+light_len}-",
                "chain='Z'")
    cmd.alter(f"{name} and chain {heavy_chain} and resi {1}-{heavy_len}",
                "chain='B'")
    cmd.alter(f"{name} and chain {heavy_chain} and resi {1+heavy_len}-",
                "chain='Z'")
else:
    if "Z" not in all_chains:
        cmd.alter(f"{name} and chain {light_chain}",
                "chain='Z'")
        cmd.alter(f"{name} and chain {heavy_chain} and resi {1}-{heavy_len}",
                "chain='B'")
        cmd.alter(f"{name} and chain Z and resi {1}-{light_len}",
                "chain='A'")
    elif "X" not in all_chains:
        cmd.alter(f"{name} and chain {light_chain}",
                "chain='X'")
        cmd.alter(f"{name} and chain {heavy_chain} and resi {1}-{heavy_len}",
                "chain='B'")
        cmd.alter(f"{name} and chain X and resi {1}-{light_len}",
                "chain='A'")
    else:
        sys.exit()
if antigen_chain != "none":
    if changed != "":
        cmd.alter(f"{name} and chain {changed}", f"chain='C'")
    else:
        cmd.alter(f"{name} and chain {antigen_chain}", f"chain='C'")
cmd.select("light", f"resi {1}-{light_len} and chain A and not alt B")
cmd.select("heavy", f"resi {1}-{heavy_len} and chain B and not alt B")
if antigen_chain != "none":
    cmd.select("antigen", "chain C and not alt B")

cmd.save(f"pdb_adjusting/{name}_heavy.pdb", "heavy")
cmd.save(f"pdb_adjusting/{name}_light.pdb", "light")
for record in SeqIO.parse(f"pdb_adjusting/{name}_heavy.pdb", "pdb-atom"):
    assert final_heavy_seq in record.seq
for record in SeqIO.parse(f"pdb_adjusting/{name}_light.pdb", "pdb-atom"):
    assert final_light_seq in record.seq
if antigen_chain != "none":
   cmd.save(f"pdb_adjusting/{name}_antigen.pdb", "antigen")

