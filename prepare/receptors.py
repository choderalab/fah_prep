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
import tempfile

import rich
import openeye
from openeye import oespruce
from openeye import oedocking
import numpy as np
from openeye import oechem
from zipfile import ZipFile

from .constants import BIOLOGICAL_SYMMETRY_HEADER, SEQRES, FRAGALYSIS_URL, MINIMUM_FRAGMENT_SIZE


# structures_path = '../Mpro_tests'
# output_basepath = '../receptors'


class PreparationConfig(NamedTuple):
    input: Path
    output: Path
    is_dimer: bool
    retain_water: Optional[bool] = False

    def __str__(self):
        msg = f"\n{str(self.input.absolute())}" \
              f"\n{str(self.output.absolute())}" \
              f"\n{str(self.is_dimer)}" \
              f"\n{str(self.retain_water)}"
        return msg


class OutputPaths(NamedTuple):
    receptor_gzipped: Path
    receptor_thiolate_gzipped: Path
    protein_pdb: Path
    protein_thiolate_pdb: Path
    ligand_pdb: Path
    ligand_sdf: Path
    ligand_mol2: Path


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


def is_diamond_structure(pdb_path: Path) -> bool:
    flag = re.search(r'-x\d+_', str(pdb_path.stem)) is not None
    return flag


def remove_from_lines(lines: List[str], string: str) -> List[str]:
    return [line for line in lines if string not in line]


def is_in_lines(lines: List[str], string: str) -> bool:
    return np.any([string in line for line in lines])


def add_prefix(lines: List[str], prefix: str) -> List[str]:
    return [line + '\n' for line in prefix.split('\n')] + lines


def lines_to_mol_graph(lines: List[str]) -> oechem.OEGraphMol:

    with tempfile.NamedTemporaryFile(delete=False, mode='wt', suffix='.pdb') as pdbfile:
        pdbfile.write(''.join(lines))
        pdbfile.close()
        complex = read_pdb_file(pdbfile.name)

    return complex


def create_molecular_graph(config: PreparationConfig) -> oechem.OEGraphMol:
    with config.input.open() as f:
        pdbfile_lines = f.readlines()

    pdbfile_lines = remove_from_lines(pdbfile_lines, 'UNK')

    if is_diamond_structure(config.input) and (not is_in_lines(pdbfile_lines, 'REMARK 350')):
        pdbfile_lines = add_prefix(pdbfile_lines, BIOLOGICAL_SYMMETRY_HEADER)

    if not config.is_dimer:
        pdbfile_lines = remove_from_lines(pdbfile_lines, 'REMARK 350')

    if not config.retain_water:
        pdbfile_lines = remove_from_lines(pdbfile_lines, 'HOH')

    pdbfile_lines = remove_from_lines(pdbfile_lines, 'LINK')

    if not is_in_lines(pdbfile_lines, 'SEQRES'):
        pdbfile_lines = add_prefix(pdbfile_lines, SEQRES)

    return lines_to_mol_graph(pdbfile_lines)


def strip_hydrogens(complex: oechem.OEGraphMol) -> oechem.OEGraphMol:
    # Strip protons from structure to allow SpruceTK to add these back
    # See: 6wnp, 6wtj, 6wtk, 6xb2, 6xqs, 6xqt, 6xqu, 6m2n
    for atom in complex.GetAtoms():
        if atom.GetAtomicNum() > 1:
            oechem.OESuppressHydrogens(atom)
    return complex


def already_prepared(outputs: OutputPaths) -> bool:
    return np.all([x.exists() for x in outputs])


def rebuild_c_terminal(complex: oechem.OEGraphMol) -> oechem.OEGraphMol:
    # Delete and rebuild C-terminal residue because Spruce causes issues with this
    # See: 6m2n 6lze
    pred = oechem.OEIsCTerminalAtom()
    for atom in complex.GetAtoms():
        if pred(atom):
            for nbor in atom.GetAtoms():
                if oechem.OEGetPDBAtomIndex(nbor) == oechem.OEPDBAtomName_O:
                    complex.DeleteAtom(nbor)
    return complex


def create_output_filenames(config: PreparationConfig) -> OutputPaths:
    prefix = config.output.joinpath(config.input.stem)
    stem = config.input.stem
    outputs = OutputPaths(
        receptor_gzipped=config.output.joinpath(f'{stem}-receptor.oeb.gz'),
        receptor_thiolate_gzipped=config.output.joinpath(f'{stem}-receptor-thiolate.oeb.gz'),
        protein_pdb=config.output.joinpath(f'{stem}-protein.pdb'),
        protein_thiolate_pdb=config.output.joinpath(f'{stem}-protein-thiolate.pdb'),
        ligand_pdb=config.output.joinpath(f'{stem}-ligand.pdb'),
        ligand_sdf=config.output.joinpath(f'{stem}-ligand.sdf'),
        ligand_mol2=config.output.joinpath(f'{stem}-ligand.mol2')
    )
    return outputs


def set_options() -> oespruce.OEMakeDesignUnitOptions:
    # Both N- and C-termini should be zwitterionic
    # Mpro cleaves its own N- and C-termini
    # See https://www.pnas.org/content/113/46/12997
    # Don't allow truncation of termini, since force fields don't have parameters for this
    # Build loops and sidechains

    opts = oespruce.OEMakeDesignUnitOptions()

    opts.GetSplitOptions().SetMinLigAtoms(MINIMUM_FRAGMENT_SIZE) # minimum fragment size (in heavy atoms)
    opts.GetPrepOptions().SetStrictProtonationMode(True)

    opts.GetPrepOptions().GetBuildOptions().SetCapNTermini(False)
    opts.GetPrepOptions().GetBuildOptions().SetCapCTermini(False)
    opts.GetPrepOptions().GetBuildOptions().SetBuildLoops(True)
    opts.GetPrepOptions().GetBuildOptions().SetBuildSidechains(True)

    opts.GetPrepOptions().GetBuildOptions().GetCapBuilderOptions().SetAllowTruncate(False)

    return opts


def prevent_flip(options: oespruce.OEMakeDesignUnitOptions, match_strings: List[str]) -> oespruce.OEMakeDesignUnitOptions:
    pred = oechem.OEAtomMatchResidue(match_strings)
    protonate_opts = options.GetPrepOptions().GetProtonateOptions()
    place_hydrogens_opts = protonate_opts.GetPlaceHydrogensOptions()
    # place_hydrogens_opts.SetBypassPredicate(pred)
    place_hydrogens_opts.SetNoFlipPredicate(pred)
    return options


def make_design_units(complex: oechem.OEGraphMol, metadata: oespruce.OEStructureMetadata,
                      options: oespruce.OEMakeDesignUnitOptions) -> List[oechem.OEDesignUnit]:
    dus = list(oespruce.OEMakeDesignUnits(complex, metadata, options))
    if len(dus) == 0:
        raise RuntimeError('No design units found.')
    return dus


# class DockingSytem:


def make_receptor(design_unit: oechem.OEDesignUnit)-> oechem.OEGraphMol:
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)
    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)

    # Create receptor and other files
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    return receptor

def prepare_receptor(config: PreparationConfig) -> None:
    output_filenames = create_output_filenames(config)

    if already_prepared(output_filenames):
        return

    complex = create_molecular_graph(config)
    complex = strip_hydrogens(complex)
    complex = rebuild_c_terminal(complex)

    # Log warnings
    errfs = oechem.oeosstream() # create a stream that writes internally to a stream
    oechem.OEThrow.SetOutputStream(errfs)
    oechem.OEThrow.Clear()
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Verbose) # capture verbose error output

    opts = set_options()
    opts = prevent_flip(opts, match_strings=["GLN:189:.*:.*:.*"])
    mdata = oespruce.OEStructureMetadata()

    design_units = make_design_units(complex, mdata, opts)
    design_unit = design_units[0]


    oedocking.OEWriteReceptorFile(receptor, str(output_filenames.receptor_gzipped))

    with oechem.oemolostream(str(output_filenames.protein_pdb)) as ofs:
        oechem.OEWriteMolecule(ofs, protein)
    with oechem.oemolostream(str(output_filenames.ligand_mol2)) as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(str(output_filenames.ligand_pdb)) as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(str(output_filenames.ligand_sdf)) as ofs:
        oechem.OEWriteMolecule(ofs, ligand)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [line for line in  output_filenames.protein_pdb.open(mode='r') if 'UNK' not in line]
    with output_filenames.protein_pdb.open(mode='wt') as outfile:
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
    oedocking.OEWriteReceptorFile(receptor, str(output_filenames.receptor_thiolate_gzipped))

    with oechem.oemolostream(str(output_filenames.protein_thiolate_pdb)) as ofs:
        oechem.OEWriteMolecule(ofs, protein)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in output_filenames.protein_thiolate_pdb.open(mode='r') if 'UNK' not in line ]
    with output_filenames.protein_thiolate_pdb.open(mode='wt') as outfile:
        outfile.write(''.join(pdbfile_lines))


def download_fragalysis_latest(structures_path: Path) -> None:
    zip_path = structures_path.joinpath('Mpro.zip')
    download_url(FRAGALYSIS_URL, zip_path)
    with ZipFile(zip_path, 'r') as zip_obj:
        zip_obj.extractall(structures_path)


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
    configs = [PreparationConfig(input=x, output=y, is_dimer=y.stem == 'dimer') for x, y in
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
            prepare_receptor(config)
