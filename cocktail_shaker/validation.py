#!/usr/bin/env python3
#
# Validation Functions for SMILES and 3D Mols
#
# -------------------------------------------


# imports
# -------
from rdkit import Chem

class RaiseMoleculeError(Exception):

    __version_error_parser__ = "1.1.0"
    __allow_update__ = False

    """

    Raise Molecule Error if for some reason we can't evaluate a SMILES, 2D, or 3D molecule.

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class MoleculeValidator(object):

    """

    Used to validate the SMILES and RDKit rendering

    """

    def __init__(self, molecules, smiles=False):

        self.molecules = molecules
        self.smiles = smiles

        if smiles:
            self.validate_smiles(self.molecules)
            self.molecules = [Chem.MolFromSmiles(molecule) for molecule in self.molecules]
            self.validate_molecule(self.molecules)
        else:
            self.validate_molecule(self.molecules)
            self.molecules = [Chem.MolToSmiles(molecule) for molecule in self.molecules]
            self.validate_smiles(self.molecules)

    def validate_smiles(self, smiles):
        """

        This method takes the smiles string and runs through the validation check of a smiles string.

        Arguments:
            self (Object): Class Cocktail
            smiles (List): List of smiles strings that need to be evaluated.
        Returns:
            N / A

        Exceptions:
            RaiseMoleculeError (Exception): MolVs Stacktrace and the smiles string that failed.

        """

        # Use MolVS to validate the smiles to make sure enumeration and r group connections are correct
        # at least in the 1D Format.
        from molvs import validate_smiles as vs

        for i in range(len(self.molecules)):
            try:
                vs(self.molecules[i])
            except RaiseMoleculeError as RME:
                print ("Not a Valid Smiles, Please check the formatting: %s" % self.original_smiles)
                print ("MolVs Stacktrace %s" % RME)

        return

    def validate_molecule(self, molecule):

        """

        This function will be used to validate molecule objects

        Arguments:
            self (Object): Class Cocktail
            molecule (RDKit Object): Molecule object we need to sanitize.
        Returns:
            molecule (RDKit Object): The RDkit Molecule molecule object

        Exceptions:
            RaiseMoleculeError (Exception): Raise the Raise Molcule Error if the molecule is not valid.

        """

        for i in range(len(self.molecules)):
            try:
                Chem.rdmolops.SanitizeMol(self.molecules[i])
            except RaiseMoleculeError as RME:
                print ("Not a valid molecule: %s" % RME)

        return molecule

class FileValidator(object):

    """

    Used to validate the SMILES and RDKit rendering

    """

    _FILE_FORMATS_ = ['alc', 'cdxml', 'cerius', 'charmm', 'cif', 'cml', 'gjf', 'gromacs', 'hyperchem', 'jme', 'mol',
                            'mol2', 'mrv', 'pdb', 'sdf3000', 'sln', 'xyz', 'sdf', 'txt']

    def __init__(self, file):

        self.file = file

    def _validate_alc(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_cdxml(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_cerius(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_charmm(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_cif(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_cml(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_gjf(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_gromacs(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_hyperchem(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_jme(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_mol(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_mol2(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_mrv(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_pdb(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_sdf3000(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_sln(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_xyz(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass

    def _validate_sdf(self):


        """

        Placeholder for alc validation of alc files.

        """
        pass