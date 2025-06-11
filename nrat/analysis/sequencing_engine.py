from ..core.structure_framework import *
from  ..data import reference_tools as ref
from ..io import file_manager as fmgr
import Bio.Align
import re
from typing import *

def alignToCanonical(chain:chain) -> Bio.Align.Alignment:
    try:
        upCode=ref.UP_Codes[chain.type()]
    except:
        raise TypeError(f"{chain.name} could not be identified")
    
    seq=fmgr.getFasta(upCode)
    aligner=Bio.Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.1
    aligner.match_score = 2
    alignments = aligner.align(chain.seq(),seq)

    best=alignments[0]
    return best

def analyzeAlignment(alignment:Bio.Align.Alignment) -> dict[str,list[tuple]]:
    '''
    Returns Gaps and Mutations in a chain alignment.
    '''
    reference,result,sample = alignment._format_unicode().split()[:3] # Slightly cursed way (but the best AFAIK) of getting the alignment values for each AA
    gaps,gap,gapStart=[],False,None
    mutations=[]

    for index in range(len(result)):
        value=result[index]
        if value == '-' and not gap:
            gap=True
            gapStart=index
        if gap and value != '-':
            gaps.append((gapStart,index))
            gap=False
        if value == '.':
            mutations.append((index+1,reference[index],sample[index]))

    return {'Gaps':gaps,'Mutations':mutations}

def fixOffset(chain:chain,alignment:Bio.Align.Alignment) -> NoReturn:
    '''
    Fixes a chain's residues to match their canonical sequence numbering scheme
    '''
    sample_index=0
    res=alignment[0]
    for x in range(len(res)):
        if res[x]!='-':
            #print(sample_index)
            chain.residues[sample_index].id=str(x+1)
            sample_index+=1


#I'm going to have this in a config later
import subprocess
blastp_path = r"C:\Program Files\blast\bin\blastp.exe"  # Change to your actual path
query_file = "query.fasta"
database = r"C:\Users\clayt\swissprot"  # Or full path to your database if not in current dir
output_file = "blast_results.txt"

def blastp(seq:str) -> tuple[str | Literal['unknown'],str | Literal['unknown']]|None:
    '''
    Protein BLAST for NR and coregulators
    '''
    query=open(query_file,'w')
    query.write('>query\n')
    query.write(seq.upper())
    query.close()
    
    cmd=[
    blastp_path,
    "-query", query_file,
    "-db", database,
    "-out", output_file,
    # sseq qseqid sseqid pident length mismatch gapopen qstart qend ssart send evalue bitscore
    "-outfmt", "6 pident stitle",     # Tabular format
    ]

    if len(seq) > 20:
        cmd+=[
            '-evalue','1e-100'
        ]
    else:
        cmd+=[
            "-evalue", "10",
            "-task","blastp-short"
        ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print("STDOUT:\n", result.stdout)
        if result.stderr:
            print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print("Error running BLASTP!")
        print("Return code:", e.returncode)
        print("Command:", e.cmd)
        print("Output:", e.output)
        print("Error output:\n", e.stderr)
    
    results=open("blast_results.txt")
    try: pident,all_names=results.readlines()[0].strip().split('\t')
    except: return None

    if float(pident)<70:
        print(f'no confident matches for given sequence {seq[:4]}...')
        return None
    
    all_names=[w.split(':')[-1].strip() for w in all_names.split(';')]
    
    filtered_names={
        'full':[],
        'short':[],
        'misc':[]
    }
    for name in all_names:
        if name.startswith('Full='):
            filtered_names['full'].append(name.split('=')[-1])
        elif name.startswith('Short='):
            filtered_names['short'].append(name.split('=')[-1])
        else:
            filtered_names['misc'].append(name)

    if len(seq)>40:
        nr_name=''
        for name in filtered_names['full']:
            if 'subfamily' in name:
                nr_name=name
        
        match = re.search(r'subfamily\s+(\d+)\s+group\s+([A-Za-z])\s+member\s+(\d+)', nr_name, re.IGNORECASE)
        if match:
            subfamily, group, member = match.groups()
            subfamily=int(subfamily)
            member=int(member)
    
        try:
            commonName=ref.NRDict[subfamily][group][member]
        except:
            return('unknown','unknown')
        
        return ('NR',commonName)
    
    else:
        try:
            coreg_name=filtered_names['short'][0]
        except:
            coreg_name=filtered_names['full'][0]
        return('coreg',coreg_name)
        