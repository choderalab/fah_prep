import filecmp
import tempfile
from pathlib import Path
from difflib import ndiff
from shutil import copy, rmtree
from typing import Optional, List

import numpy as np

from prepare.receptors import prepare_receptor, PreparationConfig


FILE_SUFFIX_TO_TEST = ['-protein.pdb','-ligand.mol2', '-ligand.sdf', '-ligand.pdb'] #,  '-protein-thiolate.pdb']
THIS_DIR = Path(__file__).parent
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


def get_n_mer(path: Path) -> str:
    if 'monomer' in path.parts:
        return 'monomer'
    else:
        return 'dimer'


def dump_bad_tests(tests: List[bool], outdir: Path, targets: List[Path]) -> None:
    n_mer = get_n_mer(targets[0])
    bad_path = data_dir.joinpath('bad', n_mer)
    for i in range(len(targets)):
        if not tests[i]:
            copy(outdir.joinpath(targets[i].name), bad_path.joinpath(targets[i].name))


def run_test(test_case: str, is_dimer: bool) -> List[bool]:
    targets = get_target_files(prefix=test_case, is_dimer=is_dimer)
    pdb_file = data_dir.joinpath(f'{test_case}.pdb')
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        prepare_receptor(PreparationConfig(input=pdb_file, output=tmp, create_dimer=is_dimer))
        tests = [files_are_same(tmp.joinpath(target.name), target) for target in targets]
        dump_bad_tests(tests, tmp, targets)
    return tests


def get_target_files(prefix: str, is_dimer: Optional[bool] = False) -> List[Path]:
    if is_dimer:
        dir = 'dimer'
    else:
        dir = 'monomer'
    return [data_dir.joinpath(dir, f"{prefix}{suffix}").absolute() for suffix in FILE_SUFFIX_TO_TEST]


def test_N_monomer() -> None:
    test_case = 'Mpro-N0029_0A_bound'
    is_dimer = False
    results = run_test(test_case, is_dimer)
    assert np.all(results)


def test_x_covalent_monomer() -> None:
    test_case = 'Mpro-x0689_0A_bound'
    is_dimer = False
    results = run_test(test_case, is_dimer)
    assert np.all(results)


def test_x_covalent_dimer() -> None:
    test_case = 'Mpro-x0689_0A_bound'
    is_dimer = True
    results = run_test(test_case, is_dimer)
    assert np.all(results)


def test_x_noncovalent_monomer() -> None:
    test_case = 'Mpro-x0072_0A_bound'
    is_dimer = False
    results = run_test(test_case, is_dimer)
    assert np.all(results)


def test_x_noncovalent_dimer() -> None:
    test_case = 'Mpro-x0072_0A_bound'
    is_dimer = True
    results = run_test(test_case, is_dimer)
    assert np.all(results)

# def test_dimer_receptor() -> None:
#     prefix = 'Mpro-N0029_0A_bound'
#     out_files = [f"{prefix}{suffix}" for suffix in FILE_SUFFIX_TO_TEST]
#     pdb_file = data_dir.joinpath(f'{prefix}.pdb')
#     with tempfile.TemporaryDirectory() as tmp:
#         tmp = Path(tmp)




