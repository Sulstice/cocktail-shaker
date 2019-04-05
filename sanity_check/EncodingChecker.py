#!/usr/bin/env python

# Imports
# ---------
from argparse import ArgumentParser

class LigandFile(object):

    """

    Experiment Converter Class, we don't have the idea required fields and optional fields just yet for experiment.

    """

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, path):

        """

        Initialize and Convert

        """

        self.path = path

        encoding = self.check_utf8()

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
            f = codecs.open(self.path, encoding='utf-8', errors='strict')
            for _ in f:
                pass
            return True
        except UnicodeDecodeError:
            return False


if __name__ == "__main__":

    parser = ArgumentParser(
        description='Pull shipping data')
    parser.add_argument("File", type=str, help="File Parser for Checker")

    args = parser.parse_args()

    LigandFile(args.File)

