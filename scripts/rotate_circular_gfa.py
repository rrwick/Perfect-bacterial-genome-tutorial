#!/usr/bin/env python3
"""
This script takes a Flye GFA file as input and it outputs the sequences (to stdout) in FASTA
format. Any circular contigs (with a single link connecting the start to the end) will be
rotated by half their length (or whatever fraction is specified by --rotate). This will serve to
move any missing/duplicated bases at the start/end of the sequence to the middle, where they can
(hopefully) be repaired by polishers.

Copyright 2022 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
import pathlib
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser(description='Rotate circular GFA sequences')
    parser.add_argument('gfa', type=pathlib.Path,
                        help='Filename of Canu assembly in GFA format')
    parser.add_argument('--rotate', type=float, default=0.5,
                        help='Rotation amount relative to the sequence length '
                             '(default: 0.5)')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    contigs, circular = load_gfa(args.gfa)
    for name, seq in contigs:
        if name in circular:
            seq, rotate_bases = rotate_seq(seq, args.rotate)
            print(f'{name}: {len(seq)} bp, circular, rotated by {rotate_bases} bp', file=sys.stderr)
        else:
            print(f'{name}: {len(seq)} bp, not circular', file=sys.stderr)
        print(f'>{name}')
        print(seq)


def rotate_seq(seq, rotate_frac):
    assert 0.0 < rotate_frac < 1.0
    start_pos = int(len(seq) * rotate_frac)
    return seq[start_pos:] + seq[:start_pos], start_pos


def load_gfa(filename):
    """
    Returns contigs as a list of (name, sequence) tuples. Also returns a set of the contig names
    which are circular (containing a circularising link and no other links).
    """
    assert str(filename).endswith('.gfa')
    contigs = []
    circular, not_circular = set(), set()
    with open(filename, 'rt') as gfa:
        for line in gfa:
            parts = line.strip().split('\t')
            if parts[0] == 'S':  # segment line
                name = parts[1]
                seq = parts[2]
                contigs.append((name, seq))
            elif parts[0] == 'L':  # link line
                seg_1, strand_1 = parts[1], parts[2]
                seg_2, strand_2 = parts[3], parts[4]
                cigar = parts[5]
                if seg_1 == seg_2 and strand_1 == strand_2 and cigar == '0M':
                    circular.add(seg_1)
                else:
                    not_circular.add(seg_1)
                    not_circular.add(seg_2)
    circular -= not_circular
    return contigs, circular


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_load_gfa():
    with tempfile.TemporaryDirectory() as tmp_dir:
        gfa_filename = pathlib.Path(tmp_dir) / 'test.gfa'
        with open(gfa_filename, 'wt') as f:
            f.write('S\tcontig1\tACGATCGACTACG\n')
            f.write('S\tcontig2\tGCCTGCCTCG\n')
            f.write('S\tcontig3\tTGGTGT\n')
            f.write('S\tcontig4\tCGCC\n')
            f.write('S\tcontig5\tTAT\n')
            f.write('S\tcontig6\tGG\n')
            f.write('S\tcontig7\tA\n')
            f.write('L\tcontig1\t+\tcontig1\t+\t0M\n')
            f.write('L\tcontig1\t-\tcontig1\t-\t0M\n')
            f.write('L\tcontig2\t+\tcontig2\t+\t0M\n')
            f.write('L\tcontig3\t-\tcontig3\t-\t0M\n')
            f.write('L\tcontig5\t+\tcontig5\t+\t0M\n')
            f.write('L\tcontig5\t-\tcontig7\t+\t0M\n')
        contigs, circular = load_gfa(gfa_filename)
        assert contigs == [('contig1', 'ACGATCGACTACG'), ('contig2', 'GCCTGCCTCG'),
                           ('contig3', 'TGGTGT'), ('contig4', 'CGCC'), ('contig5', 'TAT'),
                           ('contig6', 'GG'), ('contig7', 'A')]
        assert circular == {'contig1', 'contig2', 'contig3'}


def test_rotate_seq():
    assert rotate_seq('ACGTACGACTGG', 0.25) == 'TACGACTGGACG'
    assert rotate_seq('ACGTACGACTGG', 0.50) == 'GACTGGACGTAC'
    assert rotate_seq('ACGTACGACTGG', 0.75) == 'TGGACGTACGAC'
