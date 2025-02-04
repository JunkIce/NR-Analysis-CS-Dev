import master as PDBf
import time
import os

while True:
    query=input('CODE QUERY> ')
    print()

    smiles=PDBf.fmgr.smiles(query)
    print(smiles)
    PDBf.fmgr.pubchem_diagram(smiles,query+'.png')

    #g=PDBf.init(f'{query}.cif')
    #file=PDBf.init_struct(g)


    #file.chains['A']._chain__aminosequence=file.chains['A']._chain__aminosequence[:26]+'L'+file.chains['A']._chain__aminosequence[27:]
    #file.chains['A']._chain__aminosequence=file.chains['A']._chain__aminosequence[:45]+'L'+file.chains['A']._chain__aminosequence[46:]

    #print(file.chains['A'].align()._format_unicode().split()[1]) # slightly cursed way of extracting alignment results
    #print(PDBf.seq.analyzeAlignment(file.chains['A'].align()))


    #for x in file.chains:
        #print(file.chains[x])
    #    if len(file.chains[x].ligands)>0:
    #        print(*file.chains[x].ligands,sep='\n')
        #else:
            #print('[no associated ligands]')
        #print(f'({len(file.chains[x].waters)} Associated Water(s))')
    #print()


'''
g=PDBf.init(f'3DCT.cif')
file=PDBf.init_struct(g)
#print(file.chains['A'].align()[0])
print(file.chains['A'].residues[-1].id)
file.chains['A'].align()
print(file.chains['A'].residues[-1].id)
print(PDBf.ref.AminoacidDict[file.chains['A'].residues[-1].type])
'''
'''
files=PDBf.fmgr.csvToList(os.getcwd()+'/FXRLIST.csv')
#for x in range(len(files)):
    #files[x]=files[x]+'.cif'


mutations={}
structures=[PDBf.init(code) for code in files]
#structures=[x+'.cif' for x in structures]
structures=[PDBf.init_struct(code) for code in structures]
for file in structures:
    for id in file.chains:
        if file.chains[id].type()=='FXR':
            file.chains[id].align()
            mutations[f'{file.name}:{id}']=file.chains[id].mutations()

out=open('gaps.csv','w')
for chain in mutations:
    chainstr=chain
    for x in mutations[chain]:
        chainstr+=f',{x}'
    chainstr+='\n'
    out.write(chainstr)
out.close()
'''