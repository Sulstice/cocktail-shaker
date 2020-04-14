from peptide_builder import PeptideBuilder
from functional_group_enumerator import Cocktail
import pandas as pd

if __name__ == '__main__':


    natural_amino_acids = ["C", "CCCNC(N)=N", "CCC(N)=O", "CC(O)=O", "CS", "CCC(O)=O", "CCC(O)=O", "CCC(N)=O", "[H]",
                           "CC1=CNC=N1", "C(CC)([H])C", "CC(C)C", "CCCCN", "CCSC", "CC1=CC=CC=C1", "CO", "C(C)([H])O",
                           "CCC1=CNC2=C1C=CC=C2", "CC1=CC=C(O)C=C1", "C(C)C"]
    number = 1


    root_dataframe = pd.DataFrame(columns=["smiles", "amino_acid"])

    for amino_acid in natural_amino_acids:
        peptide_molecule = PeptideBuilder(length_of_peptide=1)
        cocktail = Cocktail(peptide_backbone=peptide_molecule,
                            ligand_library=[amino_acid],
                            enable_isomers=False)
        molecules = cocktail.shake()
        molecules = cocktail.enumerate(dimensionality='1D', enumeration_complexity='high')

        dataframe = pd.DataFrame(molecules, columns=["smiles"])

        dataframe["amino_acid"] = amino_acid
        dataframe["amino_acid_number"] = number
        number = number + 1
        root_dataframe = pd.concat([root_dataframe, dataframe])

    root_dataframe.to_hdf('data.h5', key='s', mode='w')


