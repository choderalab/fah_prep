"""
Prepare all SARS-CoV-2 Mpro structures for docking and simulation in monomer and dimer forms

This should be run from the covid-moonshot/scripts directory

"""
from typing import Optional, List, Tuple, NamedTuple
import re
import os
from pathlib import Path
import argparse
import glob
import itertools

import rich
import openeye
from openeye import oechem
from zipfile import ZipFile

from .constants import BIOLOGICAL_SYMMETRY_HEADER, SEQRES, FRAGALYSIS_URL

# structures_path = '../Mpro_tests'
# output_basepath = '../receptors'


def download_url(url, save_path, chunk_size=128):
    """
    Download file from the specified URL to the specified file path, creating base dirs if needed.
    """
    # Create directory
    import os
    base_path, filename = os.path.split(save_path)
    os.makedirs(base_path, exist_ok=True)
    # Download
    from rich.progress import track
    import requests
    r = requests.get(url, stream=True)
    with open(save_path, 'wb') as fd:
        nchunks = int(r.headers['Content-Length'])/chunk_size
        for chunk in track(r.iter_content(chunk_size=chunk_size), 'Downloading ZIP archive of Mpro structures...', total=nchunks):
            fd.write(chunk)


def read_pdb_file(pdb_file):
    #print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return (mol)


def is_diamond_structure(file_name: str) -> bool:
    flag = re.search(r'-x\d+_', file_name) is not None
    return flag


def prepare_receptor(complex_path: Path, output_path: Path,
                     dimer: Optional[bool] = False, retain_water: Optional[bool] = False) -> None:
    """
    Parameters
    ----------
    complex_path : str
        The complex PDB file to read in
    output_path : str
        Base path for output
    dimer : bool, optional, default=False
        If True, generate the dimer as the biological unit
    retain_water : bool, optional, default=False
        If True, will retain waters
    """
    # Check whether this is a diamond SARS-CoV-2 Mpro structure or not

    # complex_file_name = str(complex_file_name)
    # # output_base_path = str(output_base_path)
    # basepath = complex_file_name.parent
    # filename = complex_file_name.name
    # # basepath, filename = os.path.split(complex_file_name)
    #
    # prefix, extension = os.path.splitext(filename)
    # prefix = os.path.join(output_base_path, prefix)

    # Check if receptor already exists
    receptor_filename = output_path.joinpath(f'{complex_path.stem}-receptor.oeb.gz')
    thiolate_receptor_filename = output_path.joinpath(f'{complex_path.stem}-receptor-thiolate.oeb.gz')

    if os.path.exists(receptor_filename) and os.path.exists(thiolate_receptor_filename):
        return

    # Read in PDB file, skipping UNK atoms (left over from processing covalent ligands)
    pdbfile_lines = [line for line in open(complex_path, 'r') if 'UNK' not in line]

    # Check if biological symmetry header is present
    has_biological_symmetry_header = False
    for line in pdbfile_lines:
        if 'REMARK 350' in line:
            has_biological_symmetry_header = True
            break

    # Prepend REMARK 350 (biological symmetry) header lines for Mpro (from 5RGG) if not present
    if is_diamond_structure(complex_path) and (not has_biological_symmetry_header):
        pdbfile_lines = [ line+'\n' for line in BIOLOGICAL_SYMMETRY_HEADER.split('\n') ] + pdbfile_lines

    # If monomer is specified, drop crystal symmetry lines
    if not dimer:
        pdbfile_lines = [ line for line in pdbfile_lines if 'REMARK 350' not in line ]

    # Filter out waters
    if not retain_water:
        pdbfile_lines = [ line for line in pdbfile_lines if 'HOH' not in line ]

    # Filter out LINK records to covalent inhibitors so we can model non-covalent complex
    pdbfile_lines = [ line for line in pdbfile_lines if 'LINK' not in line ]

    # Reconstruct PDBFile contents
    pdbfile_contents = ''.join(pdbfile_lines)

    # Append SEQRES to all structures if they do not have it
    has_seqres = 'SEQRES' in pdbfile_contents
    if not has_seqres:
        #print('Adding SEQRES')
        pdbfile_contents = SEQRES + pdbfile_contents

    # Read the receptor and identify design units
    from openeye import oespruce, oechem
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(delete=False, mode='wt', suffix='.pdb') as pdbfile:
        pdbfile.write(pdbfile_contents)
        pdbfile.close()
        complex = read_pdb_file(pdbfile.name)
        # TODO: Clean up

    # Strip protons from structure to allow SpruceTK to add these back
    # See: 6wnp, 6wtj, 6wtk, 6xb2, 6xqs, 6xqt, 6xqu, 6m2n
    #print('Suppressing hydrogens')
    #print(f' Initial: {sum([1 for atom in complex.GetAtoms()])} atoms')
    for atom in complex.GetAtoms():
        if atom.GetAtomicNum() > 1:
            oechem.OESuppressHydrogens(atom)
    #print(f' Final: {sum([1 for atom in complex.GetAtoms()])} atoms')

    # Delete and rebuild C-terminal residue because Spruce causes issues with this
    # See: 6m2n 6lze
    #print('Deleting C-terminal residue O')
    pred = oechem.OEIsCTerminalAtom()
    for atom in complex.GetAtoms():
        if pred(atom):
            for nbor in atom.GetAtoms():
                if oechem.OEGetPDBAtomIndex(nbor) == oechem.OEPDBAtomName_O:
                    complex.DeleteAtom(nbor)

    #pred = oechem.OEAtomMatchResidue(["GLN:306:.*:.*:.*"])
    #for atom in complex.GetAtoms(pred):
    #    if oechem.OEGetPDBAtomIndex(atom) == oechem.OEPDBAtomName_O:
    #        print('Deleting O')
    #        complex.DeleteAtom(atom)

    #het = oespruce.OEHeterogenMetadata()
    #het.SetTitle("LIG")  # real ligand 3 letter code
    #het.SetID("CovMoonShot1234")  # in case you have corporate IDs
    #het.SetType(oespruce.OEHeterogenType_Ligand)
    #   mdata.AddHeterogenMetadata(het)

    #print('Identifying design units...')
    # Produce zero design units if we fail to protonate

    # Log warnings
    errfs = oechem.oeosstream() # create a stream that writes internally to a stream
    oechem.OEThrow.SetOutputStream(errfs)
    oechem.OEThrow.Clear()
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Verbose) # capture verbose error output

    opts = oespruce.OEMakeDesignUnitOptions()
    #print(f'ligand atoms: min {opts.GetSplitOptions().GetMinLigAtoms()}, max {opts.GetSplitOptions().GetMaxLigAtoms()}')
    opts.GetSplitOptions().SetMinLigAtoms(7) # minimum fragment size (in heavy atoms)

    mdata = oespruce.OEStructureMetadata();
    opts.GetPrepOptions().SetStrictProtonationMode(True);

    # Both N- and C-termini should be zwitterionic
    # Mpro cleaves its own N- and C-termini
    # See https://www.pnas.org/content/113/46/12997
    opts.GetPrepOptions().GetBuildOptions().SetCapNTermini(False);
    opts.GetPrepOptions().GetBuildOptions().SetCapCTermini(False);
    # Don't allow truncation of termini, since force fields don't have parameters for this
    opts.GetPrepOptions().GetBuildOptions().GetCapBuilderOptions().SetAllowTruncate(False);
    # Build loops and sidechains
    opts.GetPrepOptions().GetBuildOptions().SetBuildLoops(True);
    opts.GetPrepOptions().GetBuildOptions().SetBuildSidechains(True);

    # Don't flip Gln189
    #pred = oechem.OEAtomMatchResidue(["GLN:189: :A"])
    pred = oechem.OEAtomMatchResidue(["GLN:189:.*:.*:.*"])
    protonate_opts = opts.GetPrepOptions().GetProtonateOptions();
    place_hydrogens_opts = protonate_opts.GetPlaceHydrogensOptions()
    #place_hydrogens_opts.SetBypassPredicate(pred)
    place_hydrogens_opts.SetNoFlipPredicate(pred)
    #protonate_opts = oespruce.OEProtonateDesignUnitOptions(place_hydrogens_opts)
    #opts.GetPrepOptions().SetProtonateOptions(protonate_options);

    # Make design units
    design_units = list(oespruce.OEMakeDesignUnits(complex, mdata, opts))

    # Restore error stream
    oechem.OEThrow.SetOutputStream(oechem.oeerr)

    # Capture the warnings to a string
    warnings = errfs.str().decode("utf-8")

    if len(design_units) >= 1:
        design_unit = design_units[0]
        print('')
        print('')
        print(f'{complex_path} : SUCCESS')
        print(warnings)
    elif len(design_units) == 0:
        print('')
        print('')
        print(f'{complex_path} : FAILURE')
        print(warnings)
        msg = f'No design units found for {complex_path}\n'
        msg += warnings
        msg += '\n'
        raise Exception(msg)

    # Prepare the receptor
    #print('Preparing receptor...')
    from openeye import oedocking
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)
    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)

    # Create receptor and other files
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    oedocking.OEWriteReceptorFile(receptor, receptor_filename)

    with oechem.oemolostream(f'{prefix}-protein.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)
    with oechem.oemolostream(f'{prefix}-ligand.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.sdf') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in open(f'{prefix}-protein.pdb', 'r') if 'UNK' not in line ]
    with open(f'{prefix}-protein.pdb', 'wt') as outfile:
        outfile.write(''.join(pdbfile_lines))

    # Adjust protonation state of CYS145 to generate thiolate form
    #print('Deprotonating CYS145...') # DEBUG
    #pred = oechem.OEAtomMatchResidue(["CYS:145: :A"])
    pred = oechem.OEAtomMatchResidue(["CYS:145:.*:.*:.*"])
    place_hydrogens_opts.SetBypassPredicate(pred)
    for atom in protein.GetAtoms(pred):
        if oechem.OEGetPDBAtomIndex(atom) == oechem.OEPDBAtomName_SG:
            #print('Modifying CYS 145 SG')
            oechem.OESuppressHydrogens(atom)
            atom.SetFormalCharge(-1)
            atom.SetImplicitHCount(0)
    #print('Protonating HIS41...') # DEBUG
    #pred = oechem.OEAtomMatchResidue(["HIS:41: :A"])
    pred = oechem.OEAtomMatchResidue(["HIS:41:.*:.*:.*"])
    place_hydrogens_opts.SetBypassPredicate(pred)
    for atom in protein.GetAtoms(pred):
        if oechem.OEGetPDBAtomIndex(atom) == oechem.OEPDBAtomName_ND1:
            #print('Protonating HIS 41 ND1')
            oechem.OESuppressHydrogens(atom) # strip hydrogens from residue
            atom.SetFormalCharge(+1)
            atom.SetImplicitHCount(1)
    # Update the design unit with the modified formal charge for CYS 145 SG
    oechem.OEUpdateDesignUnit(design_unit, protein, oechem.OEDesignUnitComponents_Protein)

    # Don't flip Gln189
    #pred = oechem.OEAtomMatchResidue(["GLN:189: :A"])
    #protonate_opts = opts.GetPrepOptions().GetProtonateOptions();
    #place_hydrogens_opts = protonate_opts.GetPlaceHydrogensOptions()
    #place_hydrogens_opts.SetNoFlipPredicate(pred)

    # Adjust protonation states
    #print('Re-optimizing hydrogen positions...') # DEBUG
    #place_hydrogens_opts = oechem.OEPlaceHydrogensOptions()
    #place_hydrogens_opts.SetBypassPredicate(pred)
    #protonate_opts = oespruce.OEProtonateDesignUnitOptions(place_hydrogens_opts)
    success = oespruce.OEProtonateDesignUnit(design_unit, protonate_opts)
    design_unit.GetProtein(protein)

    # Write thiolate form of receptor
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    oedocking.OEWriteReceptorFile(receptor, thiolate_receptor_filename)

    with oechem.oemolostream(f'{prefix}-protein-thiolate.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in open(f'{prefix}-protein-thiolate.pdb', 'r') if 'UNK' not in line ]
    with open(f'{prefix}-protein-thiolate.pdb', 'wt') as outfile:
        outfile.write(''.join(pdbfile_lines))


def download_fragalysis_latest(structures_path: Path) -> None:
    zip_path = structures_path.joinpath('Mpro.zip')
    download_url(FRAGALYSIS_URL, zip_path)
    with ZipFile(zip_path, 'r') as zip_obj:
        zip_obj.extractall(structures_path)


class PreparationConfig(NamedTuple):
    input: Path
    output: Path
    is_dimer: bool

    def __str__(self):
        msg = f"\n{str(self.input.absolute())}\n{str(self.output.absolute())}\n{str(self.is_dimer)}"
        return msg


def get_structures(args: argparse.Namespace) -> List[Path]:
    if not args.structures_directory.exists():
        download_fragalysis_latest(args.structures_directory)

    file_glob = str(args.structures_directory.joinpath(args.structures_filter))
    source_pdb_files = [Path(x) for x in glob.glob(file_glob)]
    if len(source_pdb_files) == 0:
        raise RuntimeError(f'Glob path {args.structures_directory.joinpath(args.structures_filter)} '
                           f'has matched 0 files.')

    return source_pdb_files


def define_prep_configs(args: argparse.Namespace) -> List[PreparationConfig]:
    input_paths = get_structures(args)
    output_paths = [args.output_directory.joinpath(subdir) for subdir in ['monomer', 'dimer']]
    products = list(itertools.product(input_paths, output_paths))
    configs = [PreparationConfig(input=x, output=y, is_dimer=y.stem=='dimer') for x, y in
               products]
    return configs


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare bound structures for free energy calculations.')
    parser.add_argument('--structures_directory', type=Path,
                        help='Directory containing the protein/ligand bound PDB files.')
    parser.add_argument('--structures_filter', type=str,
                        default="aligned/Mpro-*_0?/Mpro-*_0?_bound.pdb",
                        help='regex filter to find PDB files in structures_directory.')
    parser.add_argument('--output_directory', type=Path,
                        help='Directory to write prepared files.')
    parser.add_argument('--dry_run', type=bool, default=False,
                        help='Only input/output destinations are printed.')
    args = parser.parse_args()

    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Quiet)

    configs = define_prep_configs(args)

    for config in configs:
        if args.dry_run:
            print(config)
        else:
            prepare_receptor(config.input, config.output, config.is_dimer)
