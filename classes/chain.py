from typing import List


class Chain:
    def __init__(self, sequence: str, annotated_sequence: List[str], type: str):
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