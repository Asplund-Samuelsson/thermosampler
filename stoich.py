#!/usr/bin/env python3

# Import modules
import sys, os
import pandas as pd
import collections
import itertools
import argparse
import re

# Import functions
from mdf import parse_equation, read_reactions

# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()

# Main code block
def main(reaction_file, outfile_name, proton_name='C00080'):
    # Load stoichiometric matrix
    S_pd = read_reactions(open(reaction_file, 'r').read(), proton_name)
    # Write stoichiometric matrix to outfile
    S_pd.to_csv(outfile_name, "\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'reactions', type=str,
        help='Load reactions.'
        )
    parser.add_argument(
        'outfile', type=str,
        help='Write stoichiometric matrix in tab-delimited format.'
        )
    parser.add_argument(
        '-H', '--proton_name', default='C00080',
        help='Name used to identify protons.'
        )
    args = parser.parse_args()
    main(args.reactions, args.outfile, args.proton_name)
