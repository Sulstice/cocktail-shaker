#!/usr/bin/env python
#
# Runs Ligand Loader
#
# ----------------------------------------------------------


# imports
# -------
import os
import numpy as np

# RDKit Imports
# -------
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import rdmolops
from rdkit.Chem.Scaffolds import MurckoScaffold
import igraph
from collections import defaultdict
from itertools import product
from rdkit.Chem import Descriptors

# Strictly Python3 TODO: Does not support Python2 Currently
# -------
from pathlib import Path

# Regex Expressions
# -------
import re
PATTERN = re.compile("R\d+")


"""

In this function we are using a couple of constants assuming you have an SDF with one molecule per SDF file. 

In this case we are using Morgan Fingerprints with a radius of 2 to capture as much neighboring atoms as possible in 
regards to the fragment of the molecule.

"""

# Constant Paths
HOME_SOURCE = Path(__file__).resolve().parents[1]
CHEM_DATA_SOURCE = "/molport_sdf_files"
SDF_PATHS = Path(str(HOME_SOURCE) + str(CHEM_DATA_SOURCE)).glob("**/*.sdf")

# Molecular Weight Distributes
MOLECULAR_WEIGHT_DISTRIBUTION = []

# RDKit Constants
FINGERPRINT_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(2)
R_GROUPS = rdRGroupDecomposition.RGroupDecompositionParameters()
R_GROUPS.removeHydrogensPostMatch = True
R_GROUPS.alignment =True
R_GROUPS.removeAllHydrogenRGroups=True



def generate_graph(molecules, threshold=0.7):
    """

    Arguments:
        molecules (Array): THe molecules to be graphed.
        threshold (Int): Similarity to the original complex molecule.

    Returns:
        graph (Object): Graph Object of molecules mapped above the threshold.

    """
    fingerprints = [FINGERPRINT_GENERATOR.GetFingerprint(molecule) for molecule in molecules]
    grid = np.zeros( (len(fingerprints), len(fingerprints)))

    graph = igraph.Graph()
    graph.add_vertices(len(molecules))

    for (molecule_1, fingerprint_1) in enumerate(fingerprints):
        for (molecule_2, fingerprint_2) in enumerate(fingerprints):
            if DataStructs.TanimotoSimilarity(fingerprint_1, fingerprint_2) >= threshold:
                graph.add_edge(molecule_1, molecule_2)

    return graph

def makebond(target, chain):
    """

    Arguments:
        target (Object):
        chain (Object):
    Returns
        new_molecule (Object): New Molecule that has the R Group.
    """

    newmol= Chem.RWMol(rdmolops.CombineMols(target, chain))
    atoms = newmol.GetAtoms()
    mapper = defaultdict(list)

    for idx, atm in enumerate(atoms):
        atom_map_num = atm.GetAtomMapNum()
        mapper[atom_map_num].append(idx)

    print (mapper)
    for idx, a_list in mapper.items():
        if len(a_list) == 2:
            print (a_list)
            atm1, atm2 = a_list
            rm_atoms = [newmol.GetAtomWithIdx(atm1),newmol.GetAtomWithIdx(atm2)]
            nbr1 = [x.GetOtherAtom(newmol.GetAtomWithIdx(atm1)) for x in newmol.GetAtomWithIdx(atm1).GetBonds()][0]
            nbr1.SetAtomMapNum(idx)
            nbr2 = [x.GetOtherAtom(newmol.GetAtomWithIdx(atm2)) for x in newmol.GetAtomWithIdx(atm2).GetBonds()][0]
            nbr2.SetAtomMapNum(idx)

    newmol.AddBond(nbr1.GetIdx(), nbr2.GetIdx(), order=Chem.rdchem.BondType.SINGLE)

    nbr1.SetAtomMapNum(0)
    nbr2.SetAtomMapNum(0)
    newmol.RemoveAtom(rm_atoms[0].GetIdx())
    newmol.RemoveAtom(rm_atoms[1].GetIdx())
    newmol = newmol.GetMol()
    return newmol


def enumerate_molecule(core, dataset, maximum_mol=10):
    """

    Arguments:
        core (Object): maximum substructure of the chemical
        r_group (Object):
        maximum_mol (Int): Maximum amount of molecules the user wishes to generate.
    Returns:

    """

    labels = list(dataset.keys())

    # Expression Check for the label
    labels = [label for label in labels if PATTERN.match(label)]

    r_groups = np.asarray([dataset[label] for label in labels])

    i, j = r_groups.shape

    combinations = [k for k in product(range(j), repeat=i)]

    resolutions = []
    for i in combinations:
        molecule = core
        for index,j in enumerate(i):
            mol = makebond(molecule, r_groups[index][j])
            MOLECULAR_WEIGHT_DISTRIBUTION.append(Descriptors.ExactMolWt(Chem.MolFromSmiles(Chem.MolToSmiles(mol))))

    raise ImportError

def get_maxium_common_substructure_smiles(molecule, maximum_common_substructure):
    """

    Arguments:
        molecule (Object): Chem mol object that comes from the img file for the scaffold
        maximum_common_substructure (Object); Scaffold For the molecule object.
    Return:
        Smiles (Img): Smiles Img.
    """
    smarts = Chem.MolFromSmarts(maximum_common_substructure.smartsString)
    match = molecule.GetSubstructMatch(smarts)
    smiles = Chem.MolFragmentToSmiles(molecule, atomsToUse=match)
    return smiles

if __name__ == '__main__':

    for path in SDF_PATHS:
        path_in_str = str(path)
        outfile = os.path.basename(os.path.normpath(path_in_str))
        sdf = path_in_str

        molecules = [m for m in Chem.SDMolSupplier(sdf)]

        for molecule in molecules:
            if molecule != None:
                AllChem.Compute2DCoords(molecule)

        if (molecules[0] == None):
            continue

        graph = generate_graph(molecules)
        blocks=graph.blocks()
        simulated_molecules_idx = sorted(list(blocks), key=lambda x: len(x), reverse=True)

        if len(simulated_molecules_idx) == 0:
            continue

        simulated_molecules = [molecules[i] for i in simulated_molecules_idx[0]]
        scaffold = [MurckoScaffold.GetScaffoldForMol(molecule) for molecule in simulated_molecules]

        maxium_common_substructure_full = rdFMCS.FindMCS(scaffold, threshold=0.7)
        maxium_common_substructure_exact = rdFMCS.FindMCS(scaffold, threshold=0.7, completeRingsOnly=True, bondCompare=rdFMCS.BondCompare.CompareOrderExact, atomCompare=rdFMCS.AtomCompare.CompareElements)
        maxium_common_substructure_any = rdFMCS.FindMCS(scaffold, threshold=0.7, ringMatchesRingOnly=True, completeRingsOnly=True, atomCompare=rdFMCS.AtomCompare.CompareAny)

        # print (Chem.MolFromSmarts(mcs1.smartsString))
        # print (Chem.MolFromSmarts(mcs2.smartsString))
        molecules_has_core = []
        core = Chem.MolFromSmarts(maxium_common_substructure_exact.smartsString)

        for molecule in molecules:
            if molecule.HasSubstructMatch(core):
                AllChem.Compute2DCoords(molecule)
                molecules_has_core.append(molecule)

        mcs_smi = get_maxium_common_substructure_smiles(molecules_has_core[0], maxium_common_substructure_exact)
        core = Chem.MolFromSmiles(mcs_smi)

        r_group = rdRGroupDecomposition.RGroupDecomposition(core, R_GROUPS)
        for molecule in molecules_has_core:
            r_group.Add(molecule)
        r_group.Process()

        # frame = pd.DataFrame(r_group.GetRGroupsAsColumns())
        #
        # frame["Smiles"] = [Chem.MolToSmiles(molecule) for molecule in molecules_has_core]
        # PandasTools.AddMoleculeColumnToFrame(frame)
        # frame = frame[["ROMol", "Smiles", "Core", "R1", "R2", "R3"]]
        # frame.head(2)

        dataset = r_group.GetRGroupsAsColumns()
        core =  Chem.RemoveHs(dataset["Core"][0])

        res = enumerate_molecule(core, dataset)
        # Draw.MolsToGridImage(res[:20])
