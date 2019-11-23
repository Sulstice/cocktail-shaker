from peptide_builders import PeptideMolecule
from functional_group_enumerator import Cocktail
from file_handler import FileWriter

if __name__ == '__main__':

    peptide_backbone = PeptideMolecule(6)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br','I', 'Cl', 'F', 'O', 'N', 'C(=O)OC']
    )

    combinations = cocktail.shake()
    FileWriter('teklklst', combinations, 'sdf')

