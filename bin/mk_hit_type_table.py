#!/usr/bin/env python

"""
Usage: mk_hit_type_table.py [options] <in_csv>

Accept a CSV file with the following columns:
    - ID
    - Hit Type

Where hit type can be: {exon, intron, opposite}

Print a table with the following columns:
    - ID
    - Exon
    - Intron
    - Opposite

Where Exon, Intron, and Opposite will have a binary representation wheather
that kind of hit type was observed for a given ID.

Examples:
    mk_hit_type_table.py data.csv
    mk_hit_type_table.py data.csv | mk_venn_diagram.py
    cut -d, -f1,6 candiate_events.csv | mk_hit_type_table.py - | mk_venn_diagram.py
"""

import argparse
from dataclasses import dataclass
import os
import sys


@dataclass
class Protein:
    id: str
    exon: bool
    intron: bool
    opposite: bool

    def __hash__(self):
        return hash(self.id)

    def to_csv(self):
        return ",".join([
            self.id,
            "1" if self.exon else "0",
            "1" if self.intron else "0",
            "1" if self.opposite else "0"
            ])


# Custom help formatter class for argparse
class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = "Usage: "
        return super(CustomHelpFormatter, self).add_usage(
            usage, actions, groups, prefix
        )

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(action.option_strings) + " " + args_string


def eprint(*args, **kwargs):
    """
    Print to stderr.
    """
    print(*args, file=sys.stderr, **kwargs)


def read_input(file_path: str):
    """
    Opens the input file or stdin if the file_path is "-".

    Parameters:
        file_path (str): The path to the input file.

    Returns:
        file_handle: A file handle to the opened file or stdin.

    Raises:
        FileNotFoundError: If the input file is not found.
    """

    if file_path == "-":
        return sys.stdin
    elif file_path.startswith("/dev/fd/"):
        fd = int(os.path.basename(file_path))
        return os.fdopen(fd, "r")
    else:
        if not os.path.isfile(file_path):
            raise FileNotFoundError

        return open(file_path, "r")


def setup_argparse() -> argparse.ArgumentParser:
    """
    Sets up the argparse instance for command-line arguments.

    Returns:
    argparse.ArgumentParser: Configured ArgumentParser instance.
    """
    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
        formatter_class=fmt,
        usage="%(prog)s [options] <in_csv>",
        add_help=False,
        description=__doc__,
    )

    # Required arguments
    required = parser.add_argument_group("Arguments")
    required.add_argument(
        "in_csv",
        metavar="in_csv",
        nargs="?",
        type=str,
        default="-",
        help="Input file path or '-' for stdin.",
    )

    # Optional arguments
    parser._optionals.title = "Options"
    # Help argument
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Display this help message and exit.",
    )

    return parser


def setup_config() -> argparse.Namespace:
    """
    Setup configuration for the script.

    Returns:
    argparse.Namespace: Parsed command-line arguments.
    """
    parser = setup_argparse()
    config = parser.parse_args()

    # Example check: if no arguments provided, print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return config


def main():
    config = setup_config()
    eprint(f"Received arguments: {config}")

    try:
        with read_input(config.in_csv) as f:
            data = f.readlines()
            eprint(f"Data read: {len(data)} lines")

    except FileNotFoundError:
        eprint(f"Error: Input file not found: {config.input}")
        sys.exit(2)


    proteins = {}
    for line in data:
        id, hit_type = line.strip().split(",")

        if id not in proteins:
            proteins[id] = Protein(id, False, False, False)

        match hit_type:
            case "exon":
                proteins[id].exon = True
            case "intron":
                proteins[id].intron = True
            case "opposite":
                proteins[id].opposite = True

    print("ID,Exon,Intron,Opposite")

    for protein in proteins.values():
        print(protein.to_csv())

if __name__ == "__main__":
    main()
