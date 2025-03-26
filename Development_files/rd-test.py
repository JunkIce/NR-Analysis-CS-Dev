from rdkit import Chem
from rdkit.Chem import Draw

import file_manager as fmgr
#from IPython.core.display import display, SVG

# Define multiple molecules using SMILES
#smiles_list = ["CCO", "CCN", "CCCl", "CCC"]  # Ethanol, Ethylamine, Chloroethane, Propane
#molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
# Generate a grid image
#img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(200, 200))

# Show the image
#img.show()
'''
ccds=[line.strip() for line in open('pparg-ligands.csv','r').readlines()]


out=open('pparg-smiles.csv','w')
for code in ccds:
    out.write(f"{code},{fmgr.smiles(code)}\n")'
'''


lines=open('pparg-smiles.csv','r').readlines()

smiles={}
for line in lines:
    line=line.strip()
    ccd,smile=line.split(',')
    smiles[ccd]=smile


molecules=[Chem.MolFromSmiles(smiles[ccd]) for ccd in smiles]

labels=[f'{ccd}: {fmgr.ccdName(ccd)}' for ccd in smiles]

# Customize drawing options
draw_options = Draw.rdMolDraw2D.MolDrawOptions()
draw_options.bondLineWidth = 1.0  # Adjust bond thickness
draw_options.additionalAtomLabelPadding = 0.1  # Adjust atom spacing
draw_options.dotsPerAngstrom = 20  # Improves resolution
draw_options.padding = 0.0

img=Draw.MolsToGridImage(molecules,molsPerRow=14,subImgSize=(400,400),useSVG=True,drawOptions=draw_options,legends=labels)
with open("molecule_grid.svg", "w") as f:
    f.write(img)