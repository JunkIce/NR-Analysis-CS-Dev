from structure_framework import *
import reference_tools as ref
import file_manager as fmgr
import Bio.Align

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

def fixOffset(chain:chain,alignment:Bio.Align.Alignment):
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