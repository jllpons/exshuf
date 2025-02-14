#!/usr/bin/env python

"""
Usage: aggr_exshuf_gffs.py [options] --id <id> --exon <exonic_gff> --intron <intronic_gff> --opposite <opp_strand_gff>

Given the ID and the 3 GFF files that the exshuf.nf pipeline produces for each gene,
this script will aggregate the data and output a CSV formatted overview of the data.

Examples:
    aggr_exshuf_gffs.py --id <id> --exon <exonic_gff> --intron <intronic_gff> --opposite <opposite_gff>
"""

import argparse
from dataclasses import dataclass
from enum import Enum
import os
import sys
from typing import (
    Dict,
    List,
    TextIO,
)


@dataclass
class GffFeature:
    seqid: str
    source: str
    _type: str
    start: int
    end: int
    score: float | str
    strand: str | str
    phase: int | str
    attributes: Dict[str, str]


class EventType(Enum):
    EXON = "exon"
    INTRON = "intron"
    OPPOSITE = "opposite"


@dataclass
class AggregatedTableRow:
    id: str
    chromosome: str
    start: int
    end: int
    strand: str
    event_type: EventType
    hmm_hit_id: str
    hmm_hit_score: float

    def to_csv(self, sep=",") -> str:
        return sep.join(
            [
                self.id,
                self.chromosome,
                str(self.start),
                str(self.end),
                self.strand,
                self.event_type.value,
                self.hmm_hit_id,
                str(self.hmm_hit_score),
            ]
        )


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
        usage="%(prog)s [options] --exon <exonic_gff> --intron <intronic_gff> --opposite <opp_strand_gff>",
        add_help=False,
        description=__doc__,
    )

    # Required arguments
    required = parser.add_argument_group("Required Arguments")
    required.add_argument(
        "--id",
        metavar="<id>",
        type=str,
        help="The gene ID.",
        required=True,
    )
    required.add_argument(
        "--exon",
        metavar="<exonic_gff>",
        type=str,
        help="Hits that intersect with the exonic regions of the gene.",
        required=True,
    )
    required.add_argument(
        "--intron",
        metavar="<intronic_gff>",
        type=str,
        help="Hits that intersect with the intronic regions of the gene.",
        required=True,
    )
    required.add_argument(
        "--opposite",
        metavar="<opp_strand_gff>",
        type=str,
        help="Hits that intersect with the gene, but are on the opposite strand.",
        required=True,
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


def parse_gff_line(line: str) -> GffFeature:
    """
    Parse a GFF line and return a GffFeature object.

    Parameters:
        line (str): A line from a GFF file.

    Returns:
        GffFeature: A GffFeature object.
    """
    fields = line.strip().split("\t")
    raw_attrs = fields[8].split(";")
    attributes = {}
    for attr in raw_attrs:
        # Ignore empty strings and ensure there's an '='
        attr = attr.strip()
        if "=" in attr:
            key_val = attr.split("=", 1)  # split only once
            if len(key_val) == 2:
                key, value = key_val
                attributes[key] = value
            else:
                # Here you can decide whether to log, raise an error, or just skip
                pass

    return GffFeature(
        seqid=fields[0],
        source=fields[1],
        _type=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=float(fields[5]) if fields[5] != "." else ".",
        strand=fields[6],
        phase=int(fields[7]) if fields[7] != "." else ".",
        attributes=attributes,
    )


def parse_gff_file(file_handle: TextIO) -> List[GffFeature]:
    """
    Parse a GFF file and return a list of GffFeature objects.

    Parameters:
        file_handle: A file handle to the GFF file.

    Returns:
        list: A list of GffFeature objects.
    """
    return [parse_gff_line(line) for line in file_handle]


def gff_to_aggregated_table_row(
    gff: GffFeature, id: str, event_type: EventType
) -> AggregatedTableRow:
    """
    Convert a GffFeature object to an AggregatedTableRow object.

    Parameters:
        gff (GffFeature): A GffFeature object.
        id (str): The gene ID.
        event_type (EventType): The type of event (exon, intron, opposite).

    Returns:
        AggregatedTableRow: An AggregatedTableRow object.
    """
    try:
        hmm_hit_id = gff.attributes["DomainName"]
    except KeyError:
        hmm_hit_id = "N/A"

    return AggregatedTableRow(
        id=id,
        chromosome=gff.seqid,
        start=gff.start,
        end=gff.end,
        strand=gff.strand,
        event_type=event_type,
        hmm_hit_id=hmm_hit_id,
        hmm_hit_score=gff.score,
    )


def main():
    config = setup_config()
    eprint(f"Received arguments: {config}")

    try:
        with read_input(config.exon) as exon_file:
            exon_gffs = parse_gff_file(exon_file)

        with read_input(config.intron) as intron_file:
            intron_gffs = parse_gff_file(intron_file)

        with read_input(config.opposite) as opposite_file:
            opposite_gffs = parse_gff_file(opposite_file)

    except FileNotFoundError:
        eprint(f"Error: Input file not found: {config.input}")
        sys.exit(2)

    aggregated_rows = []

    for gff in exon_gffs:
        aggregated_rows.append(
            gff_to_aggregated_table_row(gff, config.id, EventType.EXON)
        )

    for gff in intron_gffs:
        aggregated_rows.append(
            gff_to_aggregated_table_row(gff, config.id, EventType.INTRON)
        )

    for gff in opposite_gffs:
        aggregated_rows.append(
            gff_to_aggregated_table_row(gff, config.id, EventType.OPPOSITE)
        )

    for row in aggregated_rows:
        try:
            print(row.to_csv())

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)


if __name__ == "__main__":
    main()
