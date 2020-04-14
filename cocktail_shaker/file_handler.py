#!/usr/bin/env python
#
# File Handling parsing and writing SDF, txt files.
#
# ----------------------------------------------------------

# imports
# -------
from rdkit import Chem
from pathlib import Path
from .request_handler import CactusRequestHandler, Resolver

class FileNotSupportedError(Exception):

    __version_error_parser__ = "1.1.0"
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

    """

    __version_parser__ = "1.1.0"
    __allow_update__ = False
    _CACTUS_FILE_FORMATS = ['alc', 'cdxml', 'cerius', 'charmm', 'cif', 'cml', 'gjf', 'gromacs', 'hyperchem', 'jme', 'mol',
                     'mol2', 'mrv', 'pdb', 'sdf3000', 'sln', 'xyz']


    def __init__(self, name, molecules, option, fragmentation=None):
        self.molecules = molecules
        self.name = name
        self.option = option
        self.smiles = False
        self.fragmentation = fragmentation # Determines if they would like the SDF split into fragments.

        if self.option not in self._CACTUS_FILE_FORMATS and self.option != 'sdf' and self.option != 'txt':
            raise FileNotSupportedError(message='Cocktail Shaker does not support that file format')

        print ("Writing Files... %s%s" % (self.name, self.option))

        if self.option in self._CACTUS_FILE_FORMATS:
            method_to_call = getattr(FileWriter, 'cactus_writer')
        else:
            option_decision = self.option + "_writer"
            method_to_call = getattr(FileWriter, option_decision)

        method_to_call(self)

    def sdf_writer(self):

        """

        Arguments:
             self (Object): Parameters to write the files.

        """

        if not self.fragmentation:
            writer = Chem.SDWriter(self.name + ".sdf")
            for i in self.molecules:
                try:
                    writer.write(i)
                except:
                    continue
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

        for i in range(len(self.molecules)):
            # Our request will use smiles.
            molecule = self.molecules[i]
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

    __version_parser__ = "1.1.0"
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

        molecule = Chem.rdmolfiles.SDMolSupplier(self.file)
        molecule = Chem.MolToSmiles(molecule[0])

        return molecule

    def parse_mol(self):

        """

        Parse a Mol File

        Returns:
             molecule (RDKit Object): The mol file converted into a Mol Object

        """

        pass

    def parse_mol2(self):

        """

        parse a mol2 SYBYL File (commonly known as TRIPOS Files)

        Returns:
            molecule (RDkit Object): The mo2 file converted into a Mol Object

        """

        pass

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
