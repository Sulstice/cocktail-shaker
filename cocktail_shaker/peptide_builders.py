#!/usr/bin/env python3
#
# Peptide Molecule Builder
#
# ----------------------------------------------------------


class PeptideBackboneMolecule(object):

    """

    Building the Peptide Molecule String to pass into RDKit

    This object will be used internally for handling the peptide construction.

    """

    __version__ = '1.1.0'
    __update__ = 'False'

    def __init__(self, length_of_peptide = 1):
        self.length_of_peptide = int(length_of_peptide)
        self.peptide_base = 'NCC(NCC(O)=O)=O'
        self.build_peptide_backbone()


    def build_peptide_backbone(self):

        """

        Builds the peptide string

        Returns:
            peptide (String): Peptide molecule string of a certain length

        """

        # Return the base of the peptide for one length.
        if self.length_of_peptide == 1:
            return self.peptide_base

        peptide_n_terminus = 'NCC('
        peptide_c_terminus = 'O)=O)'

        for i in range(self.length_of_peptide):
            peptide_n_terminus += 'NCC('
            peptide_c_terminus += '=O)'


        peptide_backbone = peptide_n_terminus + peptide_c_terminus[:-1]
        print (peptide_backbone)

    def build_c_terminus(self):

        """

        Builds the C-Terminus end of the peptide molecule

        """

        pass

    def build_n_terminus(self):

        """

        Builds the N-Terminus end of the peptide molecule

        """

        pass
