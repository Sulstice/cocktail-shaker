#!/usr/bin/env python3
#
# Peptide Molecule Builder
#
# ----------------------------------------------------------

class PeptideMolecule(object):

    """

    Building the Peptide Molecule String to pass into RDKit

    This object will be used internally for handling the peptide construction.

    """

    __version__ = '1.1.0'
    __update__ = 'False'

    def __init__(self, length_of_peptide = 1, include_proline = False):
        self.length_of_peptide = int(length_of_peptide)
        self.include_proline = include_proline
        self.peptide_base = self.build_peptide_backbone()
        self.peptide_replaced = self.build_peptide_molecule_temp_replacements()

    def __repr__(self):
        return str(self.peptide_replaced)

    def build_peptide_backbone(self):

        """

        Builds the peptide string

        Returns:
            peptide (String): Peptide molecule string of a certain length
            peptide_base (String): The base peptide of 1 length amino acid.

        """

        # Return the base of the peptide for one length.
        if self.length_of_peptide == 1:
            return 'NCC(NCC(O)=O)=O'

        peptide_n_terminus = 'NCC('
        peptide_c_terminus = 'O)=O)'

        for i in range(self.length_of_peptide):
            peptide_n_terminus += 'NCC('
            peptide_c_terminus += '=O)'

        peptide_backbone = peptide_n_terminus + peptide_c_terminus[:-1]

        if self.include_proline:
            peptide_backbone = self._add_proline(peptide_backbone)

        self.peptide_base = peptide_backbone

        return peptide_backbone

    def build_peptide_molecule_temp_replacements(self):

        """

        Build a temporary functional peptide molecule string to be used for the combinations downstream.

        """
        replaced_peptides = self.peptide_base
        for i in range(1, self.length_of_peptide + 1):
            start_index_pattern = replaced_peptides.find('NCC')
            replaced_peptides = replaced_peptides[0:start_index_pattern + 1] + 'C([*:' + str(i) + '])' + replaced_peptides[start_index_pattern + 2:]

        return replaced_peptides

    def _add_proline(self, peptide_backbone):

        """

        Adds the proline on the N-Terminus

        Arguments:
            peptide_backbone (String): Peptide Backbone Chain without slots

        Returns:
            proline_peptide_backbone (String): Modified peptide_backbone with the proline addition into the n-terminus.
        """
        proline_peptide_backbone = peptide_backbone.replace("NC", "N2CCCC2", 1)

        return proline_peptide_backbone


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

