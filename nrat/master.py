# other elements of NR-AT 
from core.structure_framework import *
import io.file_manager as fmgr
import core.geometry_engine as geo
import analysis.sequencing_engine as seq
import data.reference_tools as ref

loaded_files:dict[str,structureFile]={}
groups:dict[str,list[structureFile]]={}
selected:structureFile=None
def init(target:str,name=''):
    file=fmgr.get(target)
    if not name:
        name=target.split('/')[-1]
    structure=structureFile(file)
    loaded_files[name]=structure

def group(gname:str,structs:list[structureFile]):
    if gname not in groups:
        groups[gname]=[]
        for file in structs:
            groups[gname].append(file)
    else:
        print('Group already exists.')

def group_add(gname:str,):
    pass


def old_init_struct(path:str) -> structureFile:
    return structureFile(path)
def old_init(file:str) -> str:
    return fmgr.get(file)

def old_massInit(codes:list[str]) -> list:
    paths=[]
    for code in codes:
        paths.append(old_init(code))
    return paths

#def parse(path):
#    pass