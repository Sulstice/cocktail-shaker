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
from cocktail_shaker.request_handler import CactusRequestHandler, Resolver

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
    _CACTUS_FILE_FORMATS = ['alc', 'cdxml', 'cerius', 'charmm', 'cif', 'cml', 'gjf', 'gromacs', 'hyperchem', 'jme', 'mol',
                     'mol2', 'mrv', 'pdb', 'sdf3000', 'sln', 'xyz']


    def __init__(self, name, molecules, option, fragementation=None):
        self.molecules = molecules
        self.name = name
        self.option = option
        self.fragmentation = fragementation # Determines if they would like the SDF split into fragments.

        if self.option not in self._CACTUS_FILE_FORMATS and self.option != 'sdf' and self.option != 'txt':
            raise FileNotSupportedError

        print ("Writing Files...")

        if self.option in self._CACTUS_FILE_FORMATS:
            method_to_call = getattr(FileWriter, 'cactus_writer')
        else:
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

    def cactus_writer(self):

        """

        Arguments:
            self (Object): Parameters to write the mol2 file.

        """

        """

        TODO: Incorporate TwirlyMol into a python version

        Arguments:
             molecules (List): List of molecules we need to retrieve the 3D rendering.

        Returns:
            three_d_molecules_text (List): A list of the 3D molecules representation in text.

        Exceptions:
            RequestError: Raises when the response code from cactus is not 200 (something wrong on their server side).

        """

        three_d_molecules_text = []

        bulk = False

        if len(self.molecules) > 1:
            bulk = True
        for molecule in self.molecules:
            # Our request will use smiles.
            molecule = Chem.MolToSmiles(molecule)
            url = self._construct_api_url(molecule, file_format=self.option)

            # Parse out to the request handler
            request = CactusRequestHandler(url)
            response = request.get()

            # Parse out to the resolver
            resolver = Resolver(response)
            mol2_text = resolver.cactus_mol2_resolver()
            three_d_molecules_text.append(mol2_text)

        if not self.fragmentation:
            writer = open(self.name + "." + self.option, "w")
            for i in three_d_molecules_text:
                writer.write(i)
            writer.close()
        else:
            file_count = 1
            writer = open(self.name + str(file_count) + "." + self.option, "w")
            for i in range(0, len(three_d_molecules_text)):
                if i == self.fragmentation:
                    writer.close()
                    file_count += 1
                    writer = open(self.name + str(file_count) + "." + self.option, "w")
                writer.write(three_d_molecules_text[i])
                writer.write("\n")

        writer.close()

        return three_d_molecules_text

    @staticmethod
    def _construct_api_url(molecule, file_format=None, get3d=True, **kwargs):
        """

        Return the URL for the desired API endpoint.

        Arguments:
            molecule (String): Smiles representation of the molecule we are interested in.
            get3d (Boolean): Parameter Cactus needs to handle to return the mol2 data.

        Returns:
            url (string): The full url used to request from cactus.

        """

        # imports
        # -------
        _API_BASE = 'https://cactus.nci.nih.gov/chemical/structure'

        from urllib.parse import quote, urlencode

        kwargs['format'] = file_format
        representation = 'file'
        tautomers = 'tautomers:%s' % molecule
        url = '{}/{}/{}/{}'.format(_API_BASE, quote(tautomers), representation, 'xml/')
        if get3d:
            kwargs['get3d'] = True

        url += '?%s' % urlencode(kwargs)

        return url

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
        self.molecules = self._parse_mol2_list()
        # self.generator = self._parse_file()
        # self.dataframe = self._parse_mol2()

    def _parse_mol2_list(self, delimiter="@<TRIPOS>MOLECULE"):

        """

        Borrowed Script from Derek Jones - Filed under MIT License for Open Source Software -> Acknowledgement.
        Link: https://github.com/williamdjones/deep_protein_binding/tree/10b00835024702b6d0e73092c777fed267215ca7

        Handling the read of mol2 files into molecule objects using the RDKit function MolFromMol2Block

        Arguments:
            self (Object): the mol2 file path.

        """

        def _mol2_block_fetcher(file):
            """

            Helper function to parse mol2 files into a generator object

            Arguments:
                self (Object): the mol2 file path.

            Returns:
                mol_name (Generator Object): mol2 blocks with a key mapping of mol name and mol2 block.

            """

            molname = None
            prevline = ""
            mol2 = []
            for line in file: # line will contain the molecule name followed by a newline character
                if line.startswith(delimiter) and mol2:
                    yield (molname.replace("\n", "".join(mol2)), "")
                    molname = ""
                    mol2 = []
                elif prevline.startswith(delimiter):
                    molname = line
                mol2.append(line)
                prevline = line
            if mol2:
                yield (molname, "".join(mol2))
                molname = ""

        generator = _mol2_block_fetcher(self.file)
        molecules = []
        for key, value in generator:
            print (key)
            molecules.append(Chem.MolFromMol2Block(value, sanitize=False))


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

        print (self.dataframe)
