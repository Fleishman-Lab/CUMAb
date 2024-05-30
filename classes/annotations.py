from antpack import MultiChainAnnotator
from typing import List

class AnnotatedAntibody:
    def __init__(self, sequence: str, numbering_scheme: str):
        self.sequence = sequence
        self.numbering_scheme = numbering_scheme
        self.annotated_chains = self._get_annotated_chains(numbering_scheme)
        for annotated_chain in self.annotated_chains:
            if annotated_chain.type == 'H':
                self.heavy_chain = annotated_chain
            elif annotated_chain.type == 'L':
                self.light_chain = annotated_chain

    def _get_annotated_chains(self, numbering_scheme:str):
        annotator = MultiChainAnnotator()
        annotated_chains = []
        chain_tups = annotator.analyze_seq(self.sequence, scheme=numbering_scheme)
        for chain_tup in chain_tups:
            chain = AnnotatedChain(chain_tup[0], chain_tup[1], chain_tup[3])
            annotated_chains.append(chain)

    @property
    def sequence(self):
        return self.sequence

    @property
    def heavy_chain(self):
        return self.heavy_chain
    
    @property
    def light_chain(self):
        return self.light_chain


class AnnotatedChain:
    def __init__(self, sequence:str, annotated_sequence:List[str], type:str):
        self.sequence = sequence
        self.annotated_sequence = annotated_sequence
        self.type = type

    @property
    def sequence(self):
        return self.sequence
    
    @property
    def annotated_sequence(self):
        return self.annotated_sequence
    
    @property
    def type(self):
        return self.type