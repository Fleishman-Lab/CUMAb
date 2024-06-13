from typing import List
import re 
from enum import Enum


def generate_cdr_numbers_dict():
    """
    Generate a dictionary with the Martin numbers for each CDR.
    HC CDR1: 26-35
    HC CDR2: 50-58
    HC CDR3: 94-102
    LC CDR1: 24-34
    LC CDR2: 50-56
    LC CDR3: 89-97
    """
    ranges = {
        'HC_CDR1': (26, 35),
        'HC_CDR2': (50, 58),
        'HC_CDR3': (94, 102),
        'LC_CDR1': (24, 34),
        'LC_CDR2': (50, 56),
        'LC_CDR3': (89, 97)
    }
    
    # Initialize the dictionary to store lists
    cdr_numbers_dict = {}
    
    # Create lists for each range
    for key, (start, end) in ranges.items():
        cdr_numbers_dict[key] = [str(i) for i in range(start, end + 1)]
    
    return cdr_numbers_dict

cdr_numbers_dict = generate_cdr_numbers_dict()
HC_CDR1 = cdr_numbers_dict['HC_CDR1']
HC_CDR2 = cdr_numbers_dict['HC_CDR2']
HC_CDR3 = cdr_numbers_dict['HC_CDR3']
LC_CDR1 = cdr_numbers_dict['LC_CDR1']
LC_CDR2 = cdr_numbers_dict['LC_CDR2']
LC_CDR3 = cdr_numbers_dict['LC_CDR3']

class Chain:
    def __init__(self, sequence: str, annotated_sequence: List[str], type: str):
        self.sequence: str = sequence
        self.annotated_sequence: List[str] = annotated_sequence
        self.type: str = type
        CDRs = self._get_CDRs()
        self.CDR1 = CDRs['CDR1']
        self.CDR2 = CDRs['CDR2']
        self.CDR3 = CDRs['CDR3']
        self.CDRs: List[str] = [self.CDR1, self.CDR2, self.CDR3]
        
    def _get_CDRs(self) -> List[str]:
        heavy_CDRs = {}
        light_CDRs = {}
        if self.type == "H":
            for idx, CDR in enumerate([HC_CDR1, HC_CDR2, HC_CDR3]):
                indices = find_indices(CDR, self.annotated_sequence)
                CDR_seq = "".join([self.sequence[i] for i in indices])
                heavy_CDRs[f"CDR{idx + 1}"] = CDR_seq
            return heavy_CDRs
        elif self.type == "L" or self.type == "K":
            for idx, CDR in enumerate([LC_CDR1, LC_CDR2, LC_CDR3]):
                indices = find_indices(CDR, self.annotated_sequence)
                CDR_seq = "".join([self.sequence[i] for i in indices])
                light_CDRs[f"CDR{idx + 1}"] = CDR_seq
            return light_CDRs


def find_indices(cdr: List[int], annotated_sequence: List[str]) -> dict:
    """
    Find the indices of the CDRs in the annotated sequence.
    """
    def extract_number(s) -> str:
        """
        Gets the number part from a string in the annotated sequence that could have numbers and letters.
        """
        match = re.match(r'\d+', s)
        return match.group(0)
    indices = []
    for cdr_idx in cdr:
        for idx, annotation in enumerate(annotated_sequence):
            if extract_number(annotation) == cdr_idx:
                indices.append(idx)
    return indices