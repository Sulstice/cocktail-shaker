
"""Command line interface for ``cocktail-shaker``."""

import click

__all__ = [
    'main',
]


@click.group()
def main():
    """Run the command line interface for ``cocktail-shaker``."""


if __name__ == '__main__':
    main()
