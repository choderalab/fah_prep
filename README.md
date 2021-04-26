# prepare 
Suite of functions to prepare Mpro simulations (docking, free energy) for use with folding at home. 

## receptors
Prepare receptors, protein and ligand PDB files using receptors: 
```
python prepare/receptors.py -i path/to/structures -f 'Mpro-P0047*.pdb' -o ./receptors
```
* `-i` path to structures directory
* `-f` glob filter for files
* `-o` output directory (doesn't need to exist)

**Known bugs**:

Doen't produce catalytic dyad (thiolate form)

**To-do**: 
* parallelisation
* fix thiolate bug
* covert into command line tool
* add logger
* more flexible Spruce options (currently hard-coded)
* change testu



