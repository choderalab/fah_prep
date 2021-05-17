# prepare 
Suite of functions to prepare Mpro simulations (docking, free energy) for use with folding at home. 
It has three functions: `receptors`, `ligands`, `minimize`

## Install
```
> git clone https://github.com/RobertArbon/fah_prep
> cd fah_prep
> pip install -e . 
```

## `receptors`
Prepare receptors, protein and ligand PDB files using receptors: 
```
prepare receptors -i path/to/Mpro -f 'Mpro-P0047*.pdb' -o ./receptors
```
* `-i` path to structures directory. Will download from Xchem website if not present.
  Default './MPro'. 
* `-f` glob filter for files. Should be stipulated relative to `-i`.  Default "aligned/Mpro-*_0?/Mpro-*_0?_bound.pdb" - i.e., all 
  aligned structures. 
* `-o` output directory to place the prepared files. Doesn't need to exist. Default './receptors'
* `-n` denotes dry run. In this case structures files are downloaded and all the combinations of inputs/outputs
    are printed to console. 

If the structure is a monomer, then it will create a dimer form and a monomer form. 
If the structure is a dimer, it will only produce a dimer form.

A log file is created for each structure. 

**Known bugs**:

~~Cannot control the protonation of the catalytic dyad.~~  

**To-do**: 

* Get ligand information from SDF file, not from bound protein structure. 
* parallelisation
* ~~fix thiolate bug~~
* ~~covert into command line tool~~
* more flexible Spruce options (currently hard-coded)
* change tests to reflect accurate dyads.



