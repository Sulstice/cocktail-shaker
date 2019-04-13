#!/usr/bin/env python
#
# File Handling parsing and writing SDF, txt files.
#
# ----------------------------------------------------------

# imports
# -------
from pathlib import Path

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

    TODO: Support SDF, Mol2, Mol, Smiles (TXT) File, FASTA

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

    FILE_EXTENSIONS = ['.sdf','.txt', '.text']

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

        """
        pass

    def parse_txt(self):

        """

        Parse a text file.

        """
        pass

    def parse_string(selfs):

        """

        Parse a string file.
        """

