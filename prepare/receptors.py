"""
Prepare all SARS-CoV-2 Mpro structures for docking and simulation in monomer and dimer forms

This should be run from the covid-moonshot/scripts directory

"""
import rich
import openeye

structures_path = '../Mpro_tests'
output_basepath = '../receptors'

# REMARK 350 (biological multimer symmetry operations) from 5RGG
# TODO: Remove this when REMARK 350 is restored
BIOLOGICAL_SYMMETRY_HEADER = """\
REMARK 350
REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN
REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE
REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS
REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND
REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.
REMARK 350
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC
REMARK 350 SOFTWARE USED: PISA
REMARK 350 TOTAL BURIED SURFACE AREA: 4170 ANGSTROM**2
REMARK 350 SURFACE AREA OF THE COMPLEX: 25430 ANGSTROM**2
REMARK 350 CHANGE IN SOLVENT FREE ENERGY: -3.0 KCAL/MOL
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000        0.00000
"""

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



def prepare_receptor(complex_pdb_filename, output_basepath, dimer=False, retain_water=False):
    """
    Parameters
    ----------
    complex_pdb_filename : str
        The complex PDB file to read in
    output_basepath : str
        Base path for output
    dimer : bool, optional, default=False
        If True, generate the dimer as the biological unit
    retain_water : bool, optional, default=False
        If True, will retain waters
    """
    # Check whether this is a diamond SARS-CoV-2 Mpro structure or not
    import re
    is_diamond_structure = (re.search('-x\d+_', complex_pdb_filename) is not None)
    import os
    basepath, filename = os.path.split(complex_pdb_filename)
    prefix, extension = os.path.splitext(filename)
    prefix = os.path.join(output_basepath, prefix)

    # Check if receptor already exists
    receptor_filename = f'{prefix}-receptor.oeb.gz'
    thiolate_receptor_filename = f'{prefix}-receptor-thiolate.oeb.gz'
    if os.path.exists(receptor_filename) and os.path.exists(thiolate_receptor_filename):
        return

    # Read in PDB file, skipping UNK atoms (left over from processing covalent ligands)
    pdbfile_lines = [ line for line in open(complex_pdb_filename, 'r') if 'UNK' not in line ]

    # Check if biological symmetry header is present
    has_biological_symmetry_header = False
    for line in pdbfile_lines:
        if 'REMARK 350' in line:
            has_biological_symmetry_header = True
            break

    # Prepend REMARK 350 (biological symmetry) header lines for Mpro (from 5RGG) if not present
    if is_diamond_structure and (not has_biological_symmetry_header):
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
    seqres = """\
SEQRES   1 A  306  SER GLY PHE ARG LYS MET ALA PHE PRO SER GLY LYS VAL
SEQRES   2 A  306  GLU GLY CYS MET VAL GLN VAL THR CYS GLY THR THR THR
SEQRES   3 A  306  LEU ASN GLY LEU TRP LEU ASP ASP VAL VAL TYR CYS PRO
SEQRES   4 A  306  ARG HIS VAL ILE CYS THR SER GLU ASP MET LEU ASN PRO
SEQRES   5 A  306  ASN TYR GLU ASP LEU LEU ILE ARG LYS SER ASN HIS ASN
SEQRES   6 A  306  PHE LEU VAL GLN ALA GLY ASN VAL GLN LEU ARG VAL ILE
SEQRES   7 A  306  GLY HIS SER MET GLN ASN CYS VAL LEU LYS LEU LYS VAL
SEQRES   8 A  306  ASP THR ALA ASN PRO LYS THR PRO LYS TYR LYS PHE VAL
SEQRES   9 A  306  ARG ILE GLN PRO GLY GLN THR PHE SER VAL LEU ALA CYS
SEQRES  10 A  306  TYR ASN GLY SER PRO SER GLY VAL TYR GLN CYS ALA MET
SEQRES  11 A  306  ARG PRO ASN PHE THR ILE LYS GLY SER PHE LEU ASN GLY
SEQRES  12 A  306  SER CYS GLY SER VAL GLY PHE ASN ILE ASP TYR ASP CYS
SEQRES  13 A  306  VAL SER PHE CYS TYR MET HIS HIS MET GLU LEU PRO THR
SEQRES  14 A  306  GLY VAL HIS ALA GLY THR ASP LEU GLU GLY ASN PHE TYR
SEQRES  15 A  306  GLY PRO PHE VAL ASP ARG GLN THR ALA GLN ALA ALA GLY
SEQRES  16 A  306  THR ASP THR THR ILE THR VAL ASN VAL LEU ALA TRP LEU
SEQRES  17 A  306  TYR ALA ALA VAL ILE ASN GLY ASP ARG TRP PHE LEU ASN
SEQRES  18 A  306  ARG PHE THR THR THR LEU ASN ASP PHE ASN LEU VAL ALA
SEQRES  19 A  306  MET LYS TYR ASN TYR GLU PRO LEU THR GLN ASP HIS VAL
SEQRES  20 A  306  ASP ILE LEU GLY PRO LEU SER ALA GLN THR GLY ILE ALA
SEQRES  21 A  306  VAL LEU ASP MET CYS ALA SER LEU LYS GLU LEU LEU GLN
SEQRES  22 A  306  ASN GLY MET ASN GLY ARG THR ILE LEU GLY SER ALA LEU
SEQRES  23 A  306  LEU GLU ASP GLU PHE THR PRO PHE ASP VAL VAL ARG GLN
SEQRES  24 A  306  CYS SER GLY VAL THR PHE GLN
"""
    has_seqres = 'SEQRES' in pdbfile_contents
    if not has_seqres:
        #print('Adding SEQRES')
        pdbfile_contents = seqres + pdbfile_contents

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
        print(f'{complex_pdb_filename} : SUCCESS')
        print(warnings)
    elif len(design_units) == 0:
        print('')
        print('')
        print(f'{complex_pdb_filename} : FAILURE')
        print(warnings)
        msg = f'No design units found for {complex_pdb_filename}\n'
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



if __name__ == '__main__':
    # Prep all receptors
    import glob, os

    # Be quiet
    from openeye import oechem
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Quiet)
    #oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)

    if not os.path.exists(structures_path):
        # Download ZIP file
        url = 'https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip'
        zip_path = os.path.join(structures_path, 'Mpro.zip')
        download_url(url, zip_path)
        # Unpack ZIP file
        from zipfile import ZipFile
        with ZipFile(zip_path, 'r') as zip_obj:
           zip_obj.extractall(structures_path)

    # Get list of all PDB files to prep
    source_pdb_files = glob.glob(os.path.join(structures_path, "aligned/Mpro-*_0?/Mpro-*_0?_bound.pdb"))
    #source_pdb_files = glob.glob(os.path.join(structures_path, "aligned/Mpro-6lze_0?/Mpro-*_0?_bound.pdb")) # DEBUG
    #source_pdb_files = glob.glob(os.path.join(structures_path, "aligned/Mpro-x11498_0?/Mpro-*_0?_bound.pdb")) # DEBUG

    # Create output directory
    os.makedirs(output_basepath, exist_ok=True)

    for dimer in [False, True]:
        if dimer:
            output_basepath = '../receptors/dimer'
        else:
            output_basepath = '../receptors/monomer'

        os.makedirs(output_basepath, exist_ok=True)

        def prepare_receptor_wrapper(complex_pdb_file):
            try:
                prepare_receptor(complex_pdb_file, output_basepath, dimer=dimer)
            except Exception as e:
                print(e)
                pass

        # DEBUG:
        for source_pdb_file in source_pdb_files:
           print(source_pdb_file)
           prepare_receptor_wrapper(source_pdb_file)
           print('')
        # stop
        # print(source_pdb_files)
        # Process all receptors in parallel
        # from multiprocessing import Pool
        #
        # with Pool(2) as pool:
        #     print(pool.imap_unordered(prepare_receptor_wrapper, source_pdb_files))

        # from rich.progress import Progress
        # with Pool() as pool:
        #     with Progress() as progress:
        #         task = progress.add_task('[green]Sprucing structures...', total=len(source_pdb_files))
        #         for i, _ in enumerate(pool.imap_unordered(prepare_receptor_wrapper, source_pdb_files)):
        #             progress.update(task, advance=1)