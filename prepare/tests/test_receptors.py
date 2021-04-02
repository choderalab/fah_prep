import filecmp
import tempfile
from pathlib import Path
from difflib import ndiff

import numpy as np

from prepare.receptors import prepare_receptor

FILE_SUFFIX_TO_TEST = ['-protein.pdb', '-ligand.mol2', '-ligand.pdb', '-ligand.sdf']
data_dir = Path(__file__).parent.joinpath('data')


def files_are_same(file1: Path, file2: Path) -> bool:
    if 'sdf' in file1.suffix:
        with open(file1, 'rt') as f:
            a = f.readlines()
        with open(file2, 'rt') as f:
            b = f.readlines()
        diff = list(ndiff(a, b))
        are_different = [x for x in diff if x.startswith('+') or x.startswith('-')]
        return len(are_different) == 2
    else:
        return filecmp.cmp(str(file1), str(file2))


def test_prepare_receptor() -> None:
    print('current wb', Path().absolute())

    prefix = 'Mpro-N0029_0A_bound'
    out_files = [f"{prefix}{suffix}" for suffix in FILE_SUFFIX_TO_TEST]
    pdb_file = data_dir.joinpath(f'{prefix}.pdb')
    with tempfile.TemporaryDirectory() as tmp:
        prepare_receptor(complex_pdb_filename=str(pdb_file), output_basepath=tmp)
        tmp = Path(tmp)
        tests = [files_are_same(tmp.joinpath(file), data_dir.joinpath(file)) for file in out_files]
    if not np.all(tests):
        print([out_files[i] for i in range(len(out_files)) if not tests[i]])
    assert np.all(tests)




