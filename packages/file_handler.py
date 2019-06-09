#!/usr/bin/env python
#
# File Handling parsing and writing SDF, txt files.
#
# ----------------------------------------------------------

# imports
# -------
from rdkit import Chem
from pathlib import Path
import pandas as pd

class FileNotSupportedError(Exception):

    __version_error_parser__ = 1.0
    __allow_update__ = False

    """

    Raise File Not Supported Error if we don't support the parsing. 

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class FileWriter(object):

    """

    This object is used to manage file outputs dependent on the user of the file.

    TODO: Support SDF, Mol2, Smiles (TXT) File, FASTA

    """

    __version_parser__ = 1.0
    __allow_update__ = False


    def __init__(self, name, molecules, option, fragementation=None):
        self.molecules = molecules
        self.name = name
        self.option = option
        self.fragmentation = fragementation # Determines if they would like the SDF split into fragments.

        # Avoids the continuous "if" and "else" statements.
        option_decision = self.option + "_writer"
        method_to_call = getattr(FileWriter, option_decision)
        result = method_to_call(self)

    def sdf_writer(self):

        """

        Arguments:
             self (Object): Parameters to write the files.

        """

        if not self.fragmentation:
            writer = Chem.SDWriter(self.name + ".sdf")
            for i in self.molecules:
                writer.write(i)

            writer.close()
        else:
            file_count = 1
            writer = Chem.SDWriter(self.name + str(file_count) + ".sdf")
            for i in self.molecules:
                if writer.NumMols() == self.fragmentation:
                    writer.close()
                    file_count += 1
                    writer = Chem.SDWriter(self.name + str(file_count) + ".sdf")

                writer.write(i)

            writer.close()

    def txt_writer(self):

        """

        Arguments:
             self (Object): Parameters to write the files.

        """

        if not self.fragmentation:
            writer = Chem.SmilesWriter(self.name + ".txt")
            for i in self.molecules:
                writer.write(i)

            writer.close()
        else:
            file_count = 1
            writer = Chem.SmilesWriter(self.name + str(file_count) + ".txt")
            for i in self.molecules:
                if writer.NumMols() == self.fragmentation:
                    writer.close()
                    file_count += 1
                    writer = Chem.SmilesWriter(self.name + str(file_count) + ".txt")

                writer.write(i)

            writer.close()

class FileParser(object):

    __version_parser__ = 1.0
    __allow_update__ = False

    FILE_EXTENSIONS = ['.sdf','.txt', '.text', '.mol', '.mol2']

    """

    This class function will be used to pass in files. Have to use a custom parser because we are handling chemical 
    files. 

    TODO: Support txt and SDF

    """

    def __init__(self, file):

        self.file = file
        self.check_utf8()
        file_extension = self.detect_file_extension()

        # Avoids the continuous "if" and "else" statements.
        option_decision = "parse_" + file_extension[0][1:]
        method_to_call = getattr(FileParser, option_decision)
        result = method_to_call(self)


    def check_utf8(self):

        """

        Use the codecs package to check for UTf-8 Encoding.

        Arguments:
            self (Object): Ligand File Object

        Exceptions:
            UnicodeDecodeError: Detects whether a file is utf-8 encoded using string enforce.

        """
        import codecs
        try:
            f = codecs.open(self.file, encoding='utf-8', errors='strict')
            for _ in f:
                pass
            return True
        except UnicodeDecodeError:
            return False

    def detect_file_extension(self):

        """

        Detect the File extension

        """
        suffixes = Path(self.file).suffixes

        # Handle folders
        if len(suffixes) == 0:
            print("Ligand Library Loader does not support folders in version: %s"
                  % FileParser.__version_parser__)
            raise FileNotSupportedError

        # Handle file extensions
        if len(suffixes) > 1:
            print("Ligand Library Loader does not support more than one file extension in version: %s"
                  % FileParser.__version_parser__)
            raise FileNotSupportedError

        if suffixes[0]in FileParser.FILE_EXTENSIONS:
            return suffixes
        else:
            print("Ligand Library Loader does not support more this file extension in version: %s"
                  % FileParser.__version_parser__)
            raise FileNotSupportedError

    def parse_sdf(self):

        """

        Parse an SDF file

        Returns:
            molecule (RDKit Object): The SDF file converted into a Mol Object.

        """

        molecule = Chem.SDMolSupplier(self.file)

        return molecule

    def parse_mol(self):

        """

        Parse a Mol File

        Returns:
             molecule (RDKit Object): The mol file converted into a Mol Object

        """

        molecule = Chem.MolFromMolFile(self.file)

        return molecule

    def parse_mol2(self):

        """

        parse a mol2 SYBYL File (commonly known as TRIPOS Files)

        Returns:
            molecule (RDkit Object): The mo2 file converted into a Mol Object

        """

        Mol2Parser(self.file)

    def parse_txt(self):

        """

        Parse a text file.

        """
        pass

    def parse_string(self):

        """

        Parse a strings as a potential option.

        """

        pass

class Mol2Parser(object):

    __version_parser__ = 1.0
    __allow_update__ = False

    __COLUMN_NAMES__ = ('atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id','subst_name', 'charge')
    __COLUMN_TYPES__ = (int, str, float, float, float, str, int, str, float)

    """
    
    TRIPOS Mol2 files are more complex files to handle that RDKit doesn't support. We need to handle 
    parsing of the mol2 files. 
    
    """

    def __init__(self, file):
        self.file = file
        self.generator = self._parse_file()
        self.dataframe = self._parse_mol2()
        self.molecule = self._convert_df_to_smiles()

    def _parse_file(self):

        """

        Used to the initial parsing of mol2 files and account for one or more molecules insid the file

        Arguments:
            self (Object): the mol2 file path.

        Returns:
            mol2_file (Object): Yields a generator object as a list for every mol2 file encompassed within the file.

        """

        with open(self.file, 'r') as mol2_file:
            # First item in the list will container the ID, Second item in the list will contain the row by row TRIPOS
            # molecule data
            mol2 = ['', []]
            while True:
                try:
                    row = next(mol2_file)
                    if row.startswith('@<TRIPOS>MOLECULE'):
                        if mol2[0]:
                            yield(mol2)
                        mol2 = ['', []]
                        mol2[0] = next(mol2_file).rstrip()
                        mol2[1].append(row)
                        mol2[1].append(mol2[0])
                    else:
                        mol2[1].append(row)
                except StopIteration:
                    yield (mol2)
                    return

    def _parse_mol2(self):

        """

        Load the molecule(s) into a dataframe

        Arguments:
            self (Object): the mol2 file path.

        Returns:
            df (DataFrame): the dataframe of the mol2 file first molecule.

        TODO: support multimolcule mol2 files

        """

        mol2_code, mol2_lines = next(self.generator)
        atom_section = self._atom_to_string(mol2_lines)

        df = pd.DataFrame([i.split() for i in atom_section],
                                 columns=Mol2Parser.__COLUMN_NAMES__)

        # Shape the column names
        for i in range(df.shape[1]):
            df[Mol2Parser.__COLUMN_NAMES__[i]] = df[Mol2Parser.__COLUMN_NAMES__[i]]\
                .astype(Mol2Parser.__COLUMN_TYPES__[i])

        return df

    @staticmethod
    def _atom_to_string(mol2_list):

        """

        Arguments:
            mol_list (object): list of mol2

        Returns:
            atom_section (List): a list of of mol2 strings where each row is an item in the list.

        """

        check_started = False
        first_index = 0
        last_index = 0

        for key, value  in enumerate(mol2_list):
            if value.startswith('@<TRIPOS>ATOM'):
                first_index = key + 1
                check_started = True
            elif check_started and value.startswith('@<TRIPOS>'):
                last_index = key
                break

        atom_section = mol2_list[first_index:last_index]

        return atom_section

    def _convert_df_to_smiles(self):
        """

        Convert the dataframe into a smiles string

        Arguments:
            self (Object): the mol2 file path.

        """

        from rdkit.Chem import PandasTools as PT
        SDF = PT.SaveSMILESFromFrame(self.dataframe, "test", molCol='ROMol', isomericSmiles=False)
