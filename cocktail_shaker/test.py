from peptide_builders import PeptideMolecule
from functional_group_enumerator import Cocktail

if __name__ == '__main__':

    peptide_backbone = PeptideMolecule(6)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br','I', 'Cl', 'F', 'I', 'F']
    )
    print (len(cocktail.shake()))