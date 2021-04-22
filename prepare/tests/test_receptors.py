import filecmp
import tempfile
from pathlib import Path
from difflib import ndiff
from shutil import copy

import numpy as np

from prepare.receptors import prepare_receptor, PreparationConfig


FILE_SUFFIX_TO_TEST = ['-protein.pdb','-ligand.mol2', '-ligand.sdf', '-ligand.pdb',  '-protein-thiolate.pdb']
data_dir = Path(__file__).parent.joinpath('data')


def files_are_same(file1: Path, file2: Path) -> bool:
    if ('sdf' in file1.suffix) or ('pdb' in file1.suffix):
        with open(file1, 'rt') as f:
            a = f.readlines()
        with open(file2, 'rt') as f:
            b = f.readlines()
        diff = list(ndiff(a, b))
        are_different = [x for x in diff if x.startswith('+') or x.startswith('-')]
        return len(are_different) <= 2

    else:
        return filecmp.cmp(str(file1), str(file2))


def test_prepare_receptor() -> None:
    prefix = 'Mpro-N0029_0A_bound'
    out_files = [f"{prefix}{suffix}" for suffix in FILE_SUFFIX_TO_TEST]
    pdb_file = data_dir.joinpath(f'{prefix}.pdb')
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        prepare_receptor(PreparationConfig(input=pdb_file, output=tmp, create_dimer=True))
        tests = [files_are_same(tmp.joinpath(file), data_dir.joinpath(file)) for file in out_files]

        if not np.all(tests):
            for i in range(len(out_files)):
                if not tests[i]:
                    print(out_files[i])
                    copy(tmp.joinpath(out_files[i]), Path('.').joinpath(Path(out_files[i]).name))

    assert np.all(tests)


def test_prepare_xchem_receptor() -> None:
    raise NotImplementedError


# def test_dimer_receptor() -> None:
#     prefix = 'Mpro-N0029_0A_bound'
#     out_files = [f"{prefix}{suffix}" for suffix in FILE_SUFFIX_TO_TEST]
#     pdb_file = data_dir.joinpath(f'{prefix}.pdb')
#     with tempfile.TemporaryDirectory() as tmp:
#         tmp = Path(tmp)




