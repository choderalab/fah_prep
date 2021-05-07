# prepare 
Suite of functions to prepare Mpro simulations (docking, free energy) for use with folding at home. 

## receptors
Prepare receptors, protein and ligand PDB files using receptors: 
```
python prepare/receptors.py -i Mpro/crystallographic -f 'Mpro-P0047*.pdb' -o ./receptors
```
* `-i` path to structures directory
* `-f` glob filter for files
* `-o` output directory (doesn't need to exist)

If the structure is a monomer, then it will create a dimer form and a monomer form. 
If the structure is a dimer, it will only produce a dimer form.

A log file is created for structure and dimer/monomer form. 

**Known bugs**:

Cannot control the protonation of the catalytic dyad.  

**To-do**: 

* parallelisation
* covert into command line tool
* more flexible Spruce options (currently hard-coded)
* change tests to reflect accurate dyads.



