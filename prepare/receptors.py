"""
Prepare all SARS-CoV-2 Mpro structures for docking and simulation in monomer and dimer forms

This should be run from the covid-moonshot/scripts directory

"""
from typing import Optional, List, Tuple, NamedTuple, Union
import re
from pathlib import Path
import argparse
import glob
import itertools
import tempfile

from openeye import oespruce
from openeye import oedocking
import numpy as np
from openeye import oechem
from zipfile import ZipFile

from prepare.constants import BIOLOGICAL_SYMMETRY_HEADER, SEQRES_DIMER, SEQRES_MONOMER, FRAGALYSIS_URL, MINIMUM_FRAGMENT_SIZE, CHAIN_PDB_INDEX


# structures_path = '../Mpro_tests'
# output_basepath = '../receptors'

class DockingSystem(NamedTuple):
    protein: oechem.OEGraphMol
    ligand: oechem.OEGraphMol
    receptor: oechem.OEGraphMol


class PreparationConfig(NamedTuple):
    input: Path
    output: Path
    create_dimer: bool
    retain_water: Optional[bool] = False

    def __str__(self):
        msg = f"\n{str(self.input.absolute())}" \
              f"\n{str(self.output.absolute())}" \
              f"\n{str(self.create_dimer)}" \
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

    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return mol


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


def get_chain_labels_pdb(pdb_path: Path) -> List[str]:
    pdb_lines = pdb_path.open('rt').readlines()
    return list(set([x[CHAIN_PDB_INDEX].lower() for x in pdb_lines if x.startswith('ATOM')]))


def au_is_dimeric(pdb_path: Path) -> bool:
    labels = get_chain_labels_pdb(pdb_path)
    return ('a' in labels) and ('b' in labels)


def clean_pdb(config: PreparationConfig) -> List[str]:
    with config.input.open() as f:
        pdbfile_lines = f.readlines()

    if is_diamond_structure(config.input) and (not is_in_lines(pdbfile_lines, 'REMARK 350')):
        pdbfile_lines = add_prefix(pdbfile_lines, BIOLOGICAL_SYMMETRY_HEADER)

    if not config.create_dimer:
        pdbfile_lines = remove_from_lines(pdbfile_lines, 'REMARK 350')

    if not config.retain_water:
        pdbfile_lines = remove_from_lines(pdbfile_lines, 'HOH')

    pdbfile_lines = remove_from_lines(pdbfile_lines, 'LINK')
    pdbfile_lines = remove_from_lines(pdbfile_lines, 'UNK')

    if not is_in_lines(pdbfile_lines, 'SEQRES'):
        if au_is_dimeric(config.input):
            pdbfile_lines = add_prefix(pdbfile_lines, SEQRES_DIMER)
        else:
            pdbfile_lines = add_prefix(pdbfile_lines, SEQRES_MONOMER)

    return pdbfile_lines


def strip_hydrogens(complex: oechem.OEGraphMol) -> oechem.OEGraphMol:
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

    # Prophylactic measure by JDC - not sure whether this is actually needed.
    opts = prevent_flip(opts, match_strings=["GLN:189:.*:.*:.*"])

    return opts


def prevent_flip(options: oespruce.OEMakeDesignUnitOptions, match_strings: List[str]) -> oespruce.OEMakeDesignUnitOptions:
    pred = oechem.OEAtomMatchResidue(match_strings)
    protonate_opts = options.GetPrepOptions().GetProtonateOptions()
    place_hydrogens_opts = protonate_opts.GetPlaceHydrogensOptions()
    place_hydrogens_opts.SetNoFlipPredicate(pred)
    return options


def make_design_units(complex: oechem.OEGraphMol, metadata: oespruce.OEStructureMetadata,
                      options: oespruce.OEMakeDesignUnitOptions) -> List[oechem.OEDesignUnit]:
    dus = list(oespruce.OEMakeDesignUnits(complex, metadata, options))
    if len(dus) == 0:
        raise RuntimeError('No design units found.')
    return dus


def make_docking_system(design_unit: oechem.OEDesignUnit) -> DockingSystem:
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)

    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)

    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)

    system = DockingSystem(protein=protein,
                           ligand=ligand,
                           receptor=receptor)
    return system


def write_receptor(receptor: oechem.OEGraphMol, paths: List[Path]) -> None:
    for path in paths:
        # if not path.exists():
        oedocking.OEWriteReceptorFile(oechem.OEGraphMol(receptor), str(path))


def write_protein(protein: oechem.OEGraphMol, paths: List[Path]) -> None:
    for path in paths:
        # if not path.exists():
        with oechem.oemolostream(str(path)) as ofs:
            oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(protein))

        with path.open(mode='rt') as f:
            pdbfile_lines = f.readlines()

        pdbfile_lines = remove_from_lines(pdbfile_lines, 'UNK')

        with path.open(mode='wt') as outfile:
            outfile.write(''.join(pdbfile_lines))


def write_molecular_graph(molecule: oechem.OEGraphMol, paths: List[Path]) -> None:
    for path in paths:
        # if not path.exists():
        with oechem.oemolostream(str(path)) as ofs:
            oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(molecule))


def write_docking_system(docking_system: DockingSystem, filenames: OutputPaths,
                         is_thiolate: Optional[bool] = False) -> None:
    if is_thiolate:
        receptor_path = filenames.receptor_thiolate_gzipped
        protein_path = filenames.protein_thiolate_pdb
    else:
        receptor_path = filenames.receptor_gzipped
        protein_path = filenames.protein_pdb

    write_receptor(receptor=docking_system.receptor, paths=[receptor_path])
    write_protein(protein=docking_system.protein, paths=[protein_path])
    paths = [filenames.ligand_mol2, filenames.ligand_pdb, filenames.ligand_sdf]
    write_molecular_graph(molecule=docking_system.ligand, paths=paths)


def change_protonation(molecule: oechem.OEGraphMol, options: oechem.OEPlaceHydrogensOptions,
                       match_strings: List[str], atom_name: int, formal_charge: int, implicit_h_count: int) -> \
        Tuple[oechem.OEGraphMol, oechem.OEPlaceHydrogensOptions]:

    pred = oechem.OEAtomMatchResidue(match_strings)
    for atom in molecule.GetAtoms(pred):
        if oechem.OEGetPDBAtomIndex(atom) == atom_name:
            print(f'\t ...changing protonation on {atom}')
            oechem.OESuppressHydrogens(atom)
            atom.SetImplicitHCount(implicit_h_count)
            atom.SetFormalCharge(formal_charge)
    return molecule, options


def bypass_atoms(match_strings: List[str], options: oechem.OEPlaceHydrogensOptions) -> oechem.OEPlaceHydrogensOptions:
    pred = oechem.OEAtomMatchResidue(match_strings)
    options.SetBypassPredicate(pred)
    return options


def create_dyad(state: str, docking_system: DockingSystem, design_unit: oechem.OEDesignUnit,
                    options: oespruce.OEMakeDesignUnitOptions) -> oechem.OEDesignUnit:

    protonate_opts = options.GetPrepOptions().GetProtonateOptions()
    place_h_opts = protonate_opts.GetPlaceHydrogensOptions()
    protein = docking_system.protein

    if state == 'His41(+) Cys145(-1)':
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["CYS:145:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_SG,
                                                   formal_charge=-1, implicit_h_count=0)
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["HIS:41:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_ND1,
                                                   formal_charge=+1, implicit_h_count=1)
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["HIS:41:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_NE2,
                                                   formal_charge=0, implicit_h_count=1)
    elif state == 'His41(0) Cys145(0)':
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["CYS:145:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_SG,
                                                   formal_charge=0, implicit_h_count=1)
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["HIS:41:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_ND1,
                                                   formal_charge=0, implicit_h_count=0)
        protein, place_h_opts = change_protonation(molecule=protein, options=place_h_opts,
                                                   match_strings=["HIS:41:.*:.*:.*"], atom_name=oechem.OEPDBAtomName_NE2,
                                                   formal_charge=0, implicit_h_count=1)
    else:
        ValueError("dyad_state must be one of ['His41(0) Cys145(0)', 'His41(+) Cys145(-)']")

    place_h_opts = bypass_atoms(["HIS:41:.*:.*:.*", "CYS:145:.*:.*:.*"], place_h_opts)
    oechem.OEAddExplicitHydrogens(protein)
    oechem.OEUpdateDesignUnit(design_unit, protein, oechem.OEDesignUnitComponents_Protein)
    oespruce.OEProtonateDesignUnit(design_unit, protonate_opts)
    return design_unit


def prepare_receptor(config: PreparationConfig) -> None:
    output_filenames = create_output_filenames(config)

    if au_is_dimeric(config.input) and (not config.create_dimer):
        return

    if already_prepared(output_filenames):
        return

    print('cleaning pdb...')
    pdb_lines = clean_pdb(config)
    print('creating molecular graph...')
    complex = lines_to_mol_graph(pdb_lines)

    # Some prophylatic measures by JDC - may be redundant in new OE versions
    print('stripping hydrogens...')
    complex = strip_hydrogens(complex)
    print('rebuilding c terminal...')
    complex = rebuild_c_terminal(complex)

    # Log warnings
    fname = output_filenames.protein_pdb.parent.joinpath(f"{output_filenames.protein_pdb.stem}.log")
    errfs = oechem.oeofstream(str(fname)) # create a stream that writes internally to a stream
    oechem.OEThrow.SetOutputStream(errfs)
    oechem.OEThrow.Clear()
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Debug) # capture verbose error output

    # Setup options
    print('setting options...')
    options = set_options()
    print('creating metadata...')
    metadata = oespruce.OEStructureMetadata()

    print('making base design unit...')
    design_unit = make_design_units(complex, metadata, options)[0]
    print('making base docking system...')
    docking_sytem = make_docking_system(design_unit)

    # Neutral dyad
    print('making neutral dyad...')
    print('\t updating design unit...')
    design_unit = create_dyad('His41(0) Cys145(0)', docking_sytem, design_unit, options)
    print('\t updating docking system...')
    docking_sytem = make_docking_system(design_unit)
    print('\t writing docking system...')
    write_docking_system(docking_sytem, output_filenames, is_thiolate=False)

    # Charge separated dyad
    print('making catalytic dyad...')
    print('\t updating design unit...')
    design_unit = create_dyad('His41(+) Cys145(-1)', docking_sytem, design_unit, options)
    print('\t updating docking system...')
    docking_system = make_docking_system(design_unit)
    print('\t writing docking system...')
    write_docking_system(docking_system, output_filenames, is_thiolate=True)
    errfs.close()

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
    #TODO - make logger
    [print(x) for x in source_pdb_files]
    return source_pdb_files


def define_prep_configs(args: argparse.Namespace) -> List[PreparationConfig]:
    input_paths = get_structures(args)
    output_paths = [args.output_directory.joinpath(subdir) for subdir in ['monomer', 'dimer']]
    products = list(itertools.product(input_paths, output_paths))
    configs = [PreparationConfig(input=x, output=y, create_dimer=y.stem == 'dimer') for x, y in
               products]
    return configs


def create_output_directories(configs: List[PreparationConfig]) -> None:
    for config in configs:
        if config.output.exists():
            pass
        else:
            config.output.mkdir(parents=True, exist_ok=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare bound structures for free energy calculations.')
    parser.add_argument('-i', dest='structures_directory', type=Path,
                        help='Directory containing the protein/ligand bound PDB files.')
    parser.add_argument('-f', dest='structures_filter', type=str,
                        default="aligned/Mpro-*_0?/Mpro-*_0?_bound.pdb",
                        help='glob filter to find PDB files in structures_directory.')
    parser.add_argument('-o', dest='output_directory', type=Path,
                        help='Directory to write prepared files.')
    parser.add_argument('-n', dest='dry_run', type=bool, default=False,
                        help='Only input/output destinations are printed.')
    args = parser.parse_args()

    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Quiet)

    configs = define_prep_configs(args)
    create_output_directories(configs)
    for config in configs:
        if args.dry_run:
            print(config)
        else:
            prepare_receptor(config)
