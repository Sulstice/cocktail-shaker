#!/usr/bin/env python
#
# Runs Smiles Enumerator
#
# ----------------------------------------------------------


# imports
# -------
from rdkit import Chem
from rdkit.Chem.Randomize import RandomizeMolBlock
from rdkit.Chem.ChemUtils import SDFToCSV
from rdkit.Chem import Descriptors
from progress.bar import Bar
from pathlib import Path

HOME_SOURCE = Path(__file__).resolve().parents[1]
CHEM_DATA_SOURCE = "/molport_sdf_files"
SDF_PATHS = Path(str(HOME_SOURCE) + str(CHEM_DATA_SOURCE)).glob("**/*.sdf")
MOLECULAR_WEIGHT_DISTRIBUTION = []
p = 0

class EnumerateSmiles(object):

    def __init__(self, smiles):
        self.smiles = smiles

def convert_to_sdf():
    """

    Helper Script for converting to SDF quick files.

    TODO: Will be removed into an actual functions to develop CSVs.

    """
    import csv
    import os

    for path in SDF_PATHS:
        path_in_str = str(path)
        outfile = os.path.basename(os.path.normpath(path_in_str))
        sdf = path_in_str

        suppl = Chem.SDMolSupplier(sdf)
        if suppl[0] == None:
            continue

        SDFToCSV.Convert(suppl, open("MolPort_Output.csv", "w"), keyCol=None, stopAfter=-1, includeChirality=False, smilesFrom='')

def get_molecular_weight(smiles):

    MOLECULAR_WEIGHT_DISTRIBUTION.append(Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles)))
    bar.next()

def RandomizeMol(mol):
    """

    Arguments:
        mol (RDkit Object): Mol object you would like to randomize

    Returns:
         MolBlock (RDKit Object): MolKit Object
    """
    mb = Chem.MolToMolBlock(mol)
    mb = RandomizeMolBlock(mb)
    return Chem.MolFromMolBlock(mb)

def RandomizeSmile(smile):

    """

    Arguments:
        smile (String): String of the chemical smiles
    Returns:
        smiles (String): a randomized version of that smiles string
        None

    """
    mol = RandomizeMol(Chem.MolFromSmiles(smile))
    if mol == None:
        return False
    else:
        return Chem.MolToSmiles(mol, canonical=False)

def get_mol_set(smile, tries=10000):
    """

    Arguments:
        smile (String): String of the chemical smiles
        tries (String): How many attempts to add unique smiles to a set
    Returns:
         s (Set): A set of unique smiles variations.

    """
    s = set()
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smile))
    s.add(canonical)
    for _ in range(tries):
        mol = RandomizeSmile(smile)
        if mol != None:
            # print (mol)
            s.add(mol)
    print (s)

    return s

if __name__ == "__main__":

    import csv
    bar = Bar('Processing', max=100)

    with open('MolPort_Output.csv') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            get_molecular_weight(row['SMILES'])
            with open(r"randomized_smiles.csv", 'w', newline='') as write_file:
                write=csv.writer(write_file)
                write.writerows([r] for r in list(get_mol_set(row['SMILES'])))

    with open (r'molecular_weight.csv', 'w', newline='') as write_file:
        write=csv.writer(write_file)
        write.writerows([r] for r in MOLECULAR_WEIGHT_DISTRIBUTION)

    bar.finish()
