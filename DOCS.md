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

**`structureFile(path)`**

Initialized by providing it a path to a pdb or cif file.
Initializes all child objects upon initialization.

**`chain(name, rawData, parent=None)`**

Usually initialized by a structureFile cascade.
Initializes all child objects as needed.
Derived from a Bio.PDB.Chain.Chain object.

*`chain.seq()`*
1) Initializes amino acid sequence
2) Attempts to identify self as coregulator or NR

*`chain.align()`*
Attempts to align indexing errors, stemming from gaps or etc.