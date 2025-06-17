# Biopython for parsing (Probably better than any parsing algorithm I could make)
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
import Bio.PDB.Chain
import Bio.PDB.Residue
import Bio.Align

def initStructureFile(path:str,name:str=None):
    format=path.split('.')[-1]
    name=name if name else path.split('/')[-1].split('.')[0]

    if format in ('cif','mmcif'):
        parser=MMCIFParser(QUIET=True)
    elif format in ('pdb'):
        parser=PDBParser(QUIET=True)
    else:
        raise TypeError('Invalid File Type')
    
    structure=parser.get_structure(name,path)

    for model in structure:
        for chain in model:
            

import os

initStructureFile(os.getcwd()+'\\3dct.cif')