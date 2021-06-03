# `prepare` 
Suite of functions to prepare simulations (docking, free energy) from crystal structures from [Fragalysis](https://fragalysis.diamond.ac.uk/viewer/react/landing). 
It has two functions: `receptors` (working), `omm` (not working). 

## Install

Requires `numpy`, `rich` and `openeye-toolkits` which you will have to install separately. 

```
> git clone https://github.com/RobertArbon/fah_prep
> cd fah_prep
> pip install -e . 
```

## `prepare receptors`

### Overview
Prepare receptors, protein and ligand PDB/SDF files using the [Open Eye toolkits](https://docs.eyesopen.com/toolkits/python/index.html), 
Spruce inparticular. `receptors` takes a set of crystal structures which contain ligands and prepares Spruce 'Design Units'
by building loops, adding hydrogens and determining protonation states etc. See the [Spruce](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html)
thand [OEChem/OEBio]( https://docs.eyesopen.com/toolkits/python/oechemtk/index.html#oebio-theory) documentation for more details.
The design unit, apo protein, ligand and receptor objects are then saved as output.  

### Warning
This function is currently hard-coded to work with the [Mpro](https://fragalysis.diamond.ac.uk/viewer/react/preview/target/Mpro)
data. In particular: 
1. If a data source can't be found it will download the Mpro data from a variable called `FRAGALYSIS_URL` which is set in
`prepare/constants.py` to be: https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip.
2. Mpro is a dimer, but some structures (those with an `x` in them, e.g., Mpro-x1119.pdb) have been crystallized
   as monomers.  If the structure is an `x` form , then the monomer form is created, along with a dimer form. The dimer form
   is created from the PDB `REMARK_350` in the PDB. This is added, if necessary, using the variable `REMARK_350` in
   `prepare/constants.py`. 
2. The options, particular to Mpro, for the Spruce design units are hard-coded in the function `get_options()` in `prepare/receptors.py`. 
2. The PDB remark `SEQRES` is added into the PDB if it is not present. This can be found in the variables 
   `SEQRES_MONOMER` and `SEQRES_DIMER` defined in `prepare/constants.py`.
3. Two forms of each structure are created, one with a charged catalytic dyad (with `thiolate` as a suffix) and a 
   neutral dyad between Cys-145 and His-41. 
   
### Usage

```
prepare receptors -i path/to/Mpro -f 'Mpro-P0047*.pdb' -o ./receptors
```
* `-i` path to structures directory. Will download from Xchem website if not present.
  Default `./MPro`. 
* `-f` glob filter for files. Should be stipulated relative to `-i`.  Default `aligned/Mpro-*_0?/Mpro-*_0?_bound.pdb` - i.e., all 
  aligned structures. 
* `-o` output directory to place the prepared files. Doesn't need to exist. Default `./receptors`
* `-n` denotes dry run. In this case structures files are downloaded and all the combinations of inputs/outputs
    are printed to console. 

* If the structure is a monomer, then it will create a dimer and a monomer form. 
* If the structure is a dimer, it will only produce a dimer form.
* A log file is created for each structure which captures all the Open Eye output. 
* No alignment is done. 
* All solvent molecules e.g., DMSO and water, are currently removed. There is the option for retaining the 
crystallographic waters, but this is currently hard-coded to be off.
* When alternate structures are available, only the first (as defined by Spruce) is used. 

### Known bugs/problems:
* **There are no unit tests**. 
* With dimers, the chain containing the ligand site is unpredicatable. e.g., if there is a ligand in both 
chain A and chain B then it is not clear which chain will contain the receptor site in the design unit.
* The design units created by Spruce appear to relabel chains (e.g., chain A -> chain B & chain B -> chain A). 

### To do: 

* Implement better logger. Currently uses the OpenEye logger which isn't great.   
* Use `multiprocessing` to perform calculations. Potential issue here as naive implementation with current  

## `prepare omm`
WIP. Performs OpenMM minimization. The code for this is in `prepare/dynamics.py` but is not currently integrated into
`prepare`. 




