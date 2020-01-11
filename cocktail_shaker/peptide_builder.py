#!/usr/bin/env python3
#
# Peptide Molecule Builder
#
# ----------------------------------------------------------

class CircularPeptideError(Exception):

    __version_error_parser__ = "1.1.0"
    __allow_update__ = False

    """

    Raise File Not Supported Error if we don't support the parsing. 

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class PeptideBuilder(object):

    """

    Building the Peptide Molecule String to pass into RDKit

    This object will be used internally for handling the peptide construction.

    """

    __version__ = '1.1.0'
    __update__ = False

    # Decrease RAM time.
    __slots__ = ['length_of_peptide', 'include_proline', 'circular', 'peptide_base', 'peptide_replaced']

    def __init__(self, length_of_peptide = 1, include_proline = False, circular = False):
        self.length_of_peptide = int(length_of_peptide)
        self.include_proline = include_proline
        self.circular = circular

        if self.circular:
            if self.length_of_peptide <= 2:
                print ("Circular peptides require a minimum of 3 amino acids in length")
            self.peptide_base = self.build_circular_peptide_backbone()
            self.peptide_replaced = self.build_peptide_molecule_temp_replacements()
        else:
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

    def build_circular_peptide_backbone(self):

        """

        Builds the circular peptide backbone string

        Returns:
            peptide (String): Peptide molecule string of a certain length
            peptide_base (String): The base peptide of 1 length amino acid.

        """

        # Return the base of the peptide for one length.
        if self.length_of_peptide == 3:
            return 'O=C1CNC(=O)CNC(=O)CN1'

        circular_peptide_n_terminus = 'CN1'
        circular_peptide_c_terminus = 'O=C1CNC(=O)'

        pattern = ''

        for i in range(self.length_of_peptide - 2):
            pattern += 'CNC(=O)'

        peptide_backbone = circular_peptide_c_terminus + pattern + circular_peptide_n_terminus

        self.peptide_base = peptide_backbone

        return peptide_backbone

    def build_peptide_molecule_temp_replacements(self):

        """

        Build a temporary functional peptide molecule string to be used for the combinations downstream.

        """

        replaced_peptides = self.peptide_base

        if self.circular:
            for i in range(1, self.length_of_peptide + 1):

                # Special case for circular peptides
                # Last carbon that closes the ring is pattern = CN1
                if i == self.length_of_peptide:
                    start_index_pattern = replaced_peptides.find('CN1')
                    replaced_peptides = replaced_peptides[0:start_index_pattern + 1] + '([*:' + str(i) + '])N' + replaced_peptides[start_index_pattern + 2:]
                else:
                    start_index_pattern = replaced_peptides.find('CNC')
                    replaced_peptides = replaced_peptides[0:start_index_pattern + 1] + '([*:' + str(i) + '])N' + replaced_peptides[start_index_pattern + 2:]

        else:
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

