from antpack import MultiChainAnnotator
from typing import List
from classes.chain import Chain

class Antibody:
    def __init__(self, sequence: str, numbering_scheme: str):
        self.sequence: str = sequence
        self.numbering_scheme: str = numbering_scheme
        self.chains: List[Chain] = self._get_chains()
        for chain in self.chains:
            if chain.type == "H":
                self.heavy_chain: Chain = chain
            elif chain.type == "L" or chain.type == "K":
                self.light_chain: Chain = chain

    def _get_chains(self) -> List[Chain]:
        annotator = MultiChainAnnotator(scheme=self.numbering_scheme)
        chains = []
        chain_tups = annotator.analyze_seq(self.sequence)
        for chain_tup in chain_tups:
            chain = Chain(chain_tup[0], chain_tup[1], chain_tup[3])
            chains.append(chain)
        return chains