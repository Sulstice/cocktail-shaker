# -*- coding: utf-8 -*-

"""Command line interface for ``cocktail-shaker``.

Why does this file exist, and why not put this in ``__main__``? You might
be tempted to import things from ``__main__`` later, but that will cause
problems--the code will get executed twice:

- When you run ``python3 -m cocktail_shaker`` python will execute
  ``__main__.py`` as a script. That means there won't be any
  ``cocktail_shaker.__main__`` in ``sys.modules``.
- When you import __main__ it will get executed again (as a module) because
  there's no ``cocktail_shaker.__main__`` in ``sys.modules``.

.. seealso:: http://click.pocoo.org/7/setuptools/#setuptools-integration
"""

import click

from .file_handler import FileWriter
from .functional_group_enumerator import Cocktail
from .peptide_builder import PeptideMolecule

__all__ = [
    'main',
]


@click.command()
@click.argument('length', type=int)
@click.argument('ligand-library', nargs=-1)
@click.option('--enable-isomers', is_flag=True)
@click.option('--include-amino-acids', is_flag=True)
@click.option('-o', '--output')
def main(
    length: int,
    ligand_library,
    enable_isomers: bool,
    include_amino_acids: bool,
    output,
):
    """Run the command line interface for ``cocktail-shaker``."""
    peptide_backbone = PeptideMolecule(length)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library=ligand_library,
        enable_isomers=enable_isomers,
        include_amino_acids=include_amino_acids,
    )
    molecules = cocktail.shake()

    if output:
        name, extension = output.split('.')
        FileWriter(name=name, molecules=molecules, option=extension)
    else:
        for molecule in molecules:
            click.echo(molecule)


if __name__ == '__main__':
    main()
