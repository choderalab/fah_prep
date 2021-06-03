"""
prepare: a suite of functions to prepare simulations from Fragalysis data for use with Folding@Home.
"""
import sys
from setuptools import setup, find_packages
# import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

# #########################
VERSION = '0.0.1'
ISRELEASED = True
__version__ = VERSION
# #########################

CLASSIFIERS = """\
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT 
Programming Language :: Python
Development Status :: beta 
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Unix
Operating System :: MacOS
Programming Language :: Python :: 3.8
"""

setup(
    # Self-descriptive entries which should always be present
    name='prepare',
    author='Robert Arbon',
    author_email='robert.arbon@gmail.com',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=VERSION,
    license='MIT',
    classifiers=CLASSIFIERS.splitlines(),
    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    entry_points={
        'console_scripts': [
            'prepare = prepare.main:main',
        ],
    },
    #Additional entries you may want simply uncomment the lines you want and fill in the data
    url='https://github.com/choderalab/',  # Website
    install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    platforms=['Linux',
               'Mac OS-X',
               'Unix'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.8",          # Python version restrictions

    #Manual control if final package is compressible or not, set False to prevent the .egg from being made
    zip_safe=False,
)
