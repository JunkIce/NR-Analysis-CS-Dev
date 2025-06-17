# Biopython for parsing (Probably better than any parsing algorithm I could make)
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
import Bio.PDB.Chain
import Bio.PDB.Residue
import Bio.Align
from typing import *

# Other things
import periodictable as ptable
import numpy as np
from typing import Iterable

import reference_tools as ref
import file_manager as fmgr
#import interaction_finder as inx



class structureFile:
    '''
    **Object for Protein Structure Files**

    Stores `chain` objects under the `structureFile.chains` dictionary.\n
    This object initializes all other structure file sub-objects from a given structure file passed through the `path` argument.\n
    Also attempts to set coregulator pairs, assuming conventional alternating host-coreg or coreg-host ordering of proteins.
    '''
    def __init__(self, path:str) -> NoReturn:
        self.filePath=path
        self.title=None
        self.pdbEntryID=None
        self.species=None
        self.chains:dict[str,chain]={}
        self.hbondCandidates=[]
        self.name=''.join(path.split('/')[-1].split('.')[:-1])
    
        acceptedModes=['cif','mmcif','pdbx','pdb']
        fileExt=path.split('.')[len(path.split('.'))-1]
        def readPDB():
            parser=PDBParser(QUIET=True)
            structure=parser.get_structure(self.filePath.split('/')[-1],self.filePath)
            for model in structure:
                for bio_chain in model:
                    chainIndex=bio_chain.id
                    self.chains[chainIndex] = chain(bio_chain.id,bio_chain,self.name)
        def readCif():
            parser=MMCIFParser(QUIET=True)
            structure=parser.get_structure(self.filePath.split('/')[-1],self.filePath)
            for model in structure:
                for bio_chain in model:
                    chainIndex=bio_chain.id
                    self.chains[chainIndex] = chain(bio_chain.id,bio_chain,self.name)

        if fileExt not in acceptedModes:
            raise ValueError(f'"{path}" is not a supported coordinate file type. Supported formats: {str(acceptedModes)}')
        else:
            modeIndex=acceptedModes.index(fileExt)
            if modeIndex <=2:
                mode='cif'
                readCif()
            elif modeIndex == 3:
                mode='legacy'
                readPDB()

        cache = None
        coregPairs={}
        for id in self.chains:
            type = self.chains[id].family()
            
            if cache != None:
                if self.chains[cache].family() != type:
                    coregPairs[cache] = id
                    coregPairs[id] = cache
                elif self.chains[cache].family() == type:
                    cache = id
            elif cache == None:
                cache == id
        self.coregPairs:dict[str,str] = coregPairs
    
    def __str__(self) -> str:
        return f'Structure File derived from \'{self.filePath}\''



class ligand():
    def __init__(self,rawData:Bio.PDB.Residue.Residue,parent:str=None) -> NoReturn:
        self.parent:str=parent
        self.name:str=rawData.resname
        self.atoms:dict[str,atom] = {}
        self.atomlist:list[atom] = []
        self.__rings:list[list[atom]] = None
        self.__full_name=None
        self.__iupac = None
        self.__smiles = None
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj
            self.atomlist.append(atom_obj)
    
    def __str__(self) -> str:
        return f'[Chain {self.parent}] Ligand: "{self.full_name() if self.full_name()!=self.name else self.iupac_name()}" (CCD: {self.name})'

    def full_name(self) -> str:
        if not self.__full_name:
            self.__full_name = fmgr.ccdName(self.name)
        return self.__full_name
    
    def iupac_name(self) -> str:
        if not self.__iupac:
            self.__iupac = fmgr.iupacName(self.name)
        return self.__iupac
    
    def smiles(self) -> str:
        if not self.__smiles:
            self.__smiles = fmgr.smiles(self.name)

    def rings(self) -> list[list[list]]:
        if not self.__rings:
            self.__rings=[]
            import geometry_engine as geo
            return geo.detect_aromatic_rings(self)

    


class atom():
    def __init__(self,serial:str,name:str,altLoc:str,resName:str,chainID:str,resSeq:str,iCode:str,x:str,y:str,z:str,occupancy:str,tempFactor:str,element:str,charge:str) -> NoReturn:
        self.serial=int(serial)
        self.name=name
        self.altLoc=altLoc
        self.resName=resName
        self.chainID=chainID
        self.resSeq=int(resSeq)
        self.iCode=iCode
        self.x=float(x)
        self.y=float(y)
        self.z=float(z)
        self.pos=[self.x,self.y,self.z]
        self.occupancy=float(occupancy)
        self.tempFactor=float(tempFactor)
        if len(element) == 2:
            element=element[0]+element[1].lower()
        self.element=element
        self.charge=charge
        #self.weight:float = ptable.mass.mass(ptable.elements.isotope(element))
    def __str__(self) -> str:
        return f'[{self.resSeq} {self.resName}]: {self.name}'
    


class residue():
    def __init__(self,rawData:Bio.PDB.Residue.Residue,parent:str=None) -> NoReturn:
        self.parent=parent
        self.id:str=str(rawData.id[1])
        self.auth_id=str(rawData.id[1])
        self.type:str=rawData.resname
        self.atoms:dict[str,atom] = {}
        self.__rings:list[list[atom]] = []
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj
        
        if self.type in ref.aa_pi_atoms:
            for loc in ref.aa_pi_atoms[self.type]:
                try: # In case of missing sidechain/part of sidechain
                    ring_atoms = [self.atoms[atom_name] for atom_name in loc]
                    self.__rings.append(ring_atoms)
                except: continue
    
    def rings(self) -> list[list[atom]]:
        return self.__rings
    def __str__(self) -> str:
        return f'[{self.parent}] Residue {self.id}: {self.type()}'

class water():
    def __init__(self,rawData:Bio.PDB.Residue.Residue,parent:str=None) -> NoReturn:
        self.parent=parent
        self.atoms:dict[str,atom]={}
        for bio_atom in rawData:
            atom_obj = atom(
            bio_atom.serial_number, bio_atom.name, bio_atom.altloc,
            rawData.resname, rawData.id[0], rawData.id[1],
            rawData.id[2], str(bio_atom.coord[0]), str(bio_atom.coord[1]),
            str(bio_atom.coord[2]), bio_atom.occupancy, bio_atom.bfactor,
            bio_atom.element, str(bio_atom.get_charge())
            )
            self.atoms[bio_atom.id] = atom_obj

class chain():
    def __init__(self,name:str,rawData:Bio.PDB.Chain.Chain,parent:str=None) -> NoReturn:
        self.parent=parent
        self.name=name
        self.residues: list[residue] = []
        self.ligands: list[ligand] = []
        self.waters: list[water] = []
        self.__aminosequence=None
        self.__alignment=None
        self.__mutations=None
        self.__gaps=None
        self.distances={}
        self.__bindingSurface=-1.0
        self.__family=None
        self.__type=None
        self.flags=[]
        self.charge_clamps=['h3cc','h4_1cc','h4_8cc','h12cc']
        self.__regPair=None
        self.__ligand_hydrogen_bonds ={}
        self.__ligand_pi_interactions = {}
        self.__got_h_bonds=False
        self.__got_pi_bonds=False
        
        for bio_residue in rawData:
            if bio_residue.resname in ref.AminoacidDict:
                self.residues.append(residue(bio_residue,self.name))
                #print(f"  Adding residue {bio_residue.resname} {bio_residue.id} to chain {self.name}")
            elif bio_residue.resname != 'HOH':
                self.ligands.append(ligand(bio_residue,self.name))
            else:
                self.waters.append(water(bio_residue,self.name))
    
    def __str__(self) -> str:
        return f'[{self.parent}] Chain {self.name}: {self.type()}'

        
    def seq(self) -> str:
        if not self.__aminosequence:
            self.__aminosequence = ''
            resBuffer=-1
            for residue in self.residues:
                currentRes=residue.id
                if currentRes != resBuffer:
                    self.__aminosequence += ref.AminoacidDict[residue.type]
                    resBuffer=currentRes

            found=False
        import sequencing_engine as seqr
        try:chain_fam,chain_type=seqr.blastp(self.__aminosequence)
        except:chain_fam,chain_type=['unknown']*2
        self.__family=chain_fam
        self.__type=chain_type

        return(self.__aminosequence)

    def family(self) -> str:
        if self.__family == None:
            self.seq()
        if 'unknown' in self.flags:
            return 'unknown'
        return self.__family

    def type(self) -> str:
        if self.__type == None:
            self.seq()
        if 'unknown' in self.flags:
            return 'unknown'
        return self.__type
    
    def align(self) -> Bio.Align.Alignment:
        if not self.__alignment:
            import sequencing_engine as seqr
            self.__alignment=seqr.alignToCanonical(self)
            seqr.fixOffset(self,self.__alignment)
        return self.__alignment

    
    def testForFam(self,test:str) -> bool:
        if self.type() != test:
            return (False)
        return(True)
    
    def gaps(self) -> list[tuple]:
        if not self.__gaps:
            import sequencing_engine as seqr
            analysis=seqr.analyzeAlignment(self.align())
            self.__gaps=analysis['Gaps']
            self.__mutations = analysis['Mutations']
        return self.__gaps
    
    def mutations(self) -> list[tuple]:
        if not self.__gaps:
            import sequencing_engine as seqr
            analysis=seqr.analyzeAlignment(self.align())
            self.__gaps=analysis['Gaps']
            self.__mutations = analysis['Mutations']
        return self.__mutations
    
    def ligand_interactions(self) -> dict[str,dict[str,dict]]:
        import interaction_finder as inx
        interactions={}
        for ligand in self.ligands:
            interactions[ligand.name]={}
            for res in self.residues:
                interaction_info=inx.all_interactions(res,ligand)
                if interaction_info:
                    interactions[ligand.name][res.id]=interaction_info
        return interactions