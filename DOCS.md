# General Documentation for this Tool

## Not meant to be a guide, but an overview

**Structure Framework**

The format of how biological objects are stored and handled is largely derived from BioPython, with some modifications to adapt that format for this use case
The higherarchy of structure objects is as follows:
```
        structureFile <-- Directly comes from a .cif or .pdb file
            |
          chain <-- one for Each Named Protein Chain in a structure
         ___|___
        |       |
    residue    ligand <-- Each residue in a protein chain and
        |___ ___|  Any ligand molecule associated with a chain
            |
           atom <-- each residue/ligand contains a list of these
```