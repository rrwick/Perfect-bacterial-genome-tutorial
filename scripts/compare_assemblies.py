#!/usr/bin/env python3
"""
This script produces a human-readable output showing the differences between two alternative
assemblies of a genome. The two assemblies must have the same number of contigs, the contigs must
be in the same order, and corresponding contigs must have the same strand and starting position.

It can be run like this to view the results directly in the terminal:
  compare_assemblies.py assembly_1.fasta assembly_2.fasta

Or you can store the results in a file like this:
  compare_assemblies.py assembly_1.fasta assembly_2.fasta > differences_1_vs_2.txt

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
import datetime
import edlib
import gzip
import mappy
import os
import re
import pathlib
import shutil
import subprocess
import sys
import textwrap


def main():
    args = parse_args()
    check_inputs(args)
    start_time = starting_message(args)
    assembly_1, assembly_2 = load_assemblies(args.assembly_1, args.assembly_2)
    align_sequences(assembly_1, assembly_2, args.padding, args.merge, args.aligner)
    finished_message(start_time)


def parse_args():
    description = 'R|Assembly comparison tool\n' + \
                  'https://github.com/rrwick/Perfect-bacterial-genome-tutorial'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    input_args = parser.add_argument_group('Inputs')
    input_args.add_argument('assembly_1', type=str,
                            help='First assembly to compare')
    input_args.add_argument('assembly_2', type=str,
                            help='Second assembly to compare')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--padding', type=int, default=15,
                              help='Bases of additional sequence to show before/after each change')
    setting_args.add_argument('--merge', type=int, default=30,
                              help='Changes this close are merged together in the output')
    setting_args.add_argument('--aligner', type=str, choices=['mappy', 'edlib'], default='mappy',
                              help='Aligner library: mappy has affine-gap, edlib is more robust')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def check_inputs(args):
    check_python_version()
    if not pathlib.Path(args.assembly_1).is_file():
        quit_with_error(f'Error: {args.assembly_1} is not a file')
    if not pathlib.Path(args.assembly_2).is_file():
        quit_with_error(f'Error: {args.assembly_2} is not a file')
    if args.padding < 0 or args.padding > 1000:
        quit_with_error('Error: the value of --padding must be >= 0 and <= 1000')


def starting_message(args):
    section_header('Assembly comparison tool')
    explanation('This script produces a human-readable file showing all of the differences '
                'between two alternative assemblies of a genome. This can help to spot '
                'troublesome regions of the genome with a large number of changes.')
    log('Assembly 1:')
    log(f'  {args.assembly_1}')
    log()
    log('Assembly 2:')
    log(f'  {args.assembly_2}')
    log()
    log('Settings:')
    log(f'  --padding {args.padding}')
    log(f'  --merge {args.merge}')
    log(f'  --aligner {args.aligner}')
    log()
    return datetime.datetime.now()


def load_assemblies(assembly_1_filename, assembly_2_filename):
    section_header('Loading assemblies')
    log(assembly_1_filename)
    assembly_1 = load_fasta(assembly_1_filename)
    for name, seq in assembly_1:
        log(f'  {name}: {len(seq):,} bp')
    log()
    log(assembly_2_filename)
    assembly_2 = load_fasta(assembly_2_filename)
    for name, seq in assembly_2:
        log(f'  {name}: {len(seq):,} bp')
    log()
    if len(assembly_1) != len(assembly_2):
        quit_with_error('Error: the assembly_1 and assembly_2 assemblies need to contain the same '
                        'number of sequences')
    for b, a in zip(assembly_1, assembly_2):
        assembly_1_name, assembly_1_seq = b
        assembly_1_seq_len = len(assembly_1_seq)
        assembly_2_name, assembly_2_seq = a
        assembly_2_seq_len = len(assembly_2_seq)
        if assembly_1_seq_len == 0 or assembly_2_seq_len == 0:
            quit_with_error('Error: zero-length sequences are not allowed')
        ratio = assembly_1_seq_len / assembly_2_seq_len
        if ratio < 0.9 or ratio > 1.11111111:
            quit_with_error(f'Error: {assembly_1_name} and {assembly_2_name} are too different in '
                            f'length - are the files in the same order?')
    return assembly_1, assembly_2


def align_sequences(assembly_1, assembly_2, padding, merge, aligner):
    section_header('Aligning sequences')
    longest_label = get_longest_label(assembly_1, assembly_2)
    for b, a in zip(assembly_1, assembly_2):
        assembly_1_name, assembly_1_seq = b
        assembly_2_name, assembly_2_seq = a
        output_differences(assembly_1_name, assembly_1_seq, assembly_2_name, assembly_2_seq,
                           padding, merge, longest_label, aligner)


def get_longest_label(assembly_1, assembly_2):
    longest_name, longest_seq = 0, 0
    for name, seq in assembly_1 + assembly_2:
        longest_name = max(longest_name, len(name))
        longest_seq = max(longest_seq, len(str(len(seq))))
    return longest_name + 2*longest_seq + 3


def output_differences(assembly_1_name, assembly_1_seq, assembly_2_name, assembly_2_seq, padding,
                       merge, longest_label, aligner):
    log(f'Aligning {assembly_1_name} to {assembly_2_name}:')
    assembly_1_aligned, assembly_2_aligned, differences, assembly_1_pos, assembly_2_pos, diff_pos, \
        expanded_cigar = get_aligned_seqs(assembly_1_seq, assembly_2_seq, aligner)
    if len(diff_pos) == 1:
        log(f'  1 difference')
    else:
        log(f'  {len(diff_pos):,} differences')
    log(f'  worst 10 bp identity:   {worst_window_identity(expanded_cigar, 10)}%')
    log(f'  worst 100 bp identity:  {worst_window_identity(expanded_cigar, 100)}%')
    log(f'  worst 1000 bp identity: {worst_window_identity(expanded_cigar, 1000)}%')
    log()

    aligned_len = len(assembly_1_aligned)
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)

    for start, end in diff_ranges:
        # Convert positions in alignment to positions in unaligned sequences:
        assembly_1_start, assembly_1_end = assembly_1_pos[start], assembly_1_pos[end]
        assembly_2_start, assembly_2_end = assembly_2_pos[start], assembly_2_pos[end]

        # Sanity check:
        assert assembly_1_aligned[start:end].replace('-', '') == \
            assembly_1_seq[assembly_1_start:assembly_1_end]
        assert assembly_2_aligned[start:end].replace('-', '') == \
            assembly_2_seq[assembly_2_start:assembly_2_end]

        # Add 1 to starts to convert from 0-based exclusive ranges to 1-based inclusive ranges.
        assembly_1_label = f'{assembly_1_name} {assembly_1_start+1}-{assembly_1_end}:'
        assembly_2_label = f'{assembly_2_name} {assembly_2_start+1}-{assembly_2_end}:'
        assert len(assembly_1_label) <= longest_label
        assert len(assembly_2_label) <= longest_label
        assembly_1_label = assembly_1_label.rjust(longest_label)
        assembly_2_label = assembly_2_label.rjust(longest_label)

        print(f'{assembly_1_label}', assembly_1_aligned[start:end], flush=True)
        print(f'{assembly_2_label}', assembly_2_aligned[start:end], flush=True)
        print(' ' * longest_label, differences[start:end], flush=True)
        print(flush=True)
    log()


def make_diff_ranges(diff_pos, padding, merge, aligned_len):
    diff_ranges = []
    last_diff_pos = None
    for p in diff_pos:
        start = max(0, p-padding)
        end = min(aligned_len, p+padding+1)
        if not last_diff_pos:  # this is the first diff
            diff_ranges.append((start, end))
        elif p - last_diff_pos <= merge:   # this diff is close to the previous diff
            prev_start = diff_ranges[-1][0]
            diff_ranges.pop()
            diff_ranges.append((prev_start, end))
        else:   # this diff is far from the previous diff
            diff_ranges.append((start, end))
        last_diff_pos = p
    return diff_ranges


def get_aligned_seqs(assembly_1_seq, assembly_2_seq, aligner):
    cigar = get_cigar(assembly_1_seq, assembly_2_seq, aligner)
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    assembly_1_aligned, assembly_2_aligned, differences = [], [], []
    assembly_1_positions, assembly_2_positions, diff_positions = [], [], []
    new_expanded_cigar = []
    for c in expanded_cigar:
        assembly_1_positions.append(i)
        assembly_2_positions.append(j)
        if c == 'M':
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            if b_1 == b_2:
                diff = ' '
                new_expanded_cigar.append('=')
            else:
                diff = '*'
                new_expanded_cigar.append('X')
                diff_positions.append(len(differences))
            i += 1
            j += 1
        elif c == '=':
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            diff = ' '
            new_expanded_cigar.append('=')
            i += 1
            j += 1
            assert b_1 == b_2
        elif c == 'X':
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            diff = '*'
            new_expanded_cigar.append('X')
            diff_positions.append(len(differences))
            i += 1
            j += 1
            assert b_1 != b_2
        elif c == 'I':
            b_1 = assembly_1_seq[i]
            b_2 = '-'
            diff = '*'
            new_expanded_cigar.append('I')
            diff_positions.append(len(differences))
            i += 1
        elif c == 'D':
            b_1 = '-'
            b_2 = assembly_2_seq[j]
            diff = '*'
            new_expanded_cigar.append('D')
            diff_positions.append(len(differences))
            j += 1
        else:
            assert False
        assembly_1_aligned.append(b_1)
        assembly_2_aligned.append(b_2)
        differences.append(diff)
    assembly_1_positions.append(i)
    assembly_2_positions.append(j)
    assert i == len(assembly_1_seq)
    assert j == len(assembly_2_seq)
    assembly_1_aligned = ''.join(assembly_1_aligned)
    assembly_2_aligned = ''.join(assembly_2_aligned)
    differences = ''.join(differences)
    for p in diff_positions:
        assert differences[p] == '*'
    assert assembly_1_aligned.replace('-', '') == assembly_1_seq
    assert assembly_2_aligned.replace('-', '') == assembly_2_seq
    new_expanded_cigar = ''.join(new_expanded_cigar)
    return assembly_1_aligned, assembly_2_aligned, differences, \
        assembly_1_positions, assembly_2_positions, diff_positions, new_expanded_cigar


def get_cigar(assembly_1_seq, assembly_2_seq, aligner):
    if assembly_1_seq == assembly_2_seq:
        return f'{len(assembly_1_seq)}='
    if aligner == 'mappy':
        return get_cigar_with_mappy(assembly_1_seq, assembly_2_seq)
    elif aligner == 'edlib':
        return get_cigar_with_edlib(assembly_1_seq, assembly_2_seq)
    else:
        assert False


def worst_window_identity(expanded_cigar, window_size):
    cigar_len = len(expanded_cigar)
    if cigar_len <= window_size:
        return 100.0 * expanded_cigar.count('=') / cigar_len
    start, end = 0, window_size
    window_match_count = expanded_cigar[start:end].count('=')
    min_match_count = window_match_count
    while end < cigar_len:
        if expanded_cigar[start] == '=':
            window_match_count -= 1
        if expanded_cigar[end] == '=':
            window_match_count += 1
        if window_match_count < min_match_count:
            min_match_count = window_match_count
        start += 1
        end += 1
    return 100.0 * min_match_count / window_size


def get_cigar_with_mappy(assembly_1_seq, assembly_2_seq):
    a = mappy.Aligner(seq=assembly_2_seq, preset='map-ont')
    for result in a.map(assembly_1_seq):
        full_length_query = (result.q_st == 0 and result.q_en == len(assembly_1_seq))
        full_length_ref = (result.r_st == 0 and result.r_en == len(assembly_2_seq))
        if full_length_query and full_length_ref:
            return result.cigar_str
    quit_with_error('Error: mappy alignment failed, try using --aligner edlib')


def get_cigar_with_edlib(assembly_1_seq, assembly_2_seq):
    result = edlib.align(assembly_1_seq, assembly_2_seq, mode='NW', task='path')
    return result['cigar']


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=M]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def finished_message(start_time):
    section_header('Finished!')
    time_to_run = datetime.datetime.now() - start_time
    log(f'Time to run: {time_to_run}')
    log()


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(fasta_filename):
    fasta_seqs = []
    with get_open_func(fasta_filename)(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if a command is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 2:  # if a command was given but nothing else
            self.print_help(file=sys.stderr)
            sys.exit(2)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default' not in help_text.lower():
                help_text += ' (default: {})'.format(action.default)
            elif 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT',
                                              'default: {}'.format(action.default))
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _split_lines(self, text, width):
        """
        Override this method to add special behaviour for help texts that start with:
          'R|' - loop text one option per line
        """
        if text.startswith('R|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2
                    line_parts = line.split(', ')
                    join = ','
                    current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith('R|'):
            return ''.join(indent + line for line in text[2:].splitlines(keepends=True))
        else:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        # determine the required width and the entry label
        help_position = min(self._action_max_length + 2,
                            self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)

        # ho nelp; start on same line and add a final newline
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0

        # short action name; start on the same line and pad two spaces
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0

        # long action name; start on the next line
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position

        # collect the pieces of the action help
        parts = [action_header]

        # if there was help for the action, add lines of help text
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))

        # or add a newline if the description doesn't end with one
        elif not action_header.endswith('\n'):
            parts.append('\n')

        # if there are any sub-actions, add their help as well
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))

        # return a single string
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


def log(message='', end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)


def section_header(text):
    log()
    time = get_timestamp()
    time_str = dim('(' + time + ')')
    header = bold_yellow_underline(text)
    print(header + ' ' + time_str, file=sys.stderr, flush=True)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def explanation(text, indent_size=4):
    text = ' ' * indent_size + text
    terminal_width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(dim(line))
    log()


def quit_with_error(text):
    terminal_width, _ = get_terminal_size_stderr()
    log()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(line)
    log()
    sys.exit()


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit('\nError: compare_assemblies.py requires Python 3.6 or later')


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_get_longest_label():
    assembly_1 = [('seq1', 'ACGT'), ('seq2', 'ACGTACGTACGTACGT')]
    assembly_2 = [('seq1_polished', 'ACGT'), ('seq2_polished', 'ACGT')]
    assert get_longest_label(assembly_1, assembly_2) == 20


def test_get_expanded_cigar():
    assert get_expanded_cigar('5=') == '====='
    assert get_expanded_cigar('3=2X4=1I6=3D3=') == '===XX====I======DDD==='


def test_make_diff_ranges():
    diff_pos = [100, 110]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 121)]

    diff_pos = [100, 120]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 131)]

    diff_pos = [100, 121]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (111, 132)]

    diff_pos = [100, 150]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (140, 161)]

    diff_pos = [2, 195]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(0, 13), (185, 200)]


def test_worst_window_identity():
    import pytest
    ex = get_expanded_cigar

    assert worst_window_identity(ex('1000=1X1000='), 10) == pytest.approx(90.0)
    assert worst_window_identity(ex('1000=1X1000='), 100) == pytest.approx(99.0)
    assert worst_window_identity(ex('1000=1X1000='), 1000) == pytest.approx(99.9)

    assert worst_window_identity(ex('1X1000='), 10) == pytest.approx(90.0)
    assert worst_window_identity(ex('1X1000='), 100) == pytest.approx(99.0)
    assert worst_window_identity(ex('1X1000='), 1000) == pytest.approx(99.9)

    assert worst_window_identity(ex('1000=1X'), 10) == pytest.approx(90.0)
    assert worst_window_identity(ex('1000=1X'), 100) == pytest.approx(99.0)
    assert worst_window_identity(ex('1000=1X'), 1000) == pytest.approx(99.9)

    assert worst_window_identity(ex('1000=1X10=1X1000='), 1000) == pytest.approx(99.8)
    assert worst_window_identity(ex('1000=1X100=1X1000='), 1000) == pytest.approx(99.8)
    assert worst_window_identity(ex('1000=1X998=1X1000='), 1000) == pytest.approx(99.8)
    assert worst_window_identity(ex('1000=1X999=1X1000='), 1000) == pytest.approx(99.9)

    assert worst_window_identity(ex('1000=5I10=5D1000='), 1000) == pytest.approx(99.0)
    assert worst_window_identity(ex('1000=5I100=5D1000='), 1000) == pytest.approx(99.0)
    assert worst_window_identity(ex('1000=5I1000=5D1000='), 1000) == pytest.approx(99.5)

    assert worst_window_identity(ex('50=1X49='), 10) == pytest.approx(90.0)
    assert worst_window_identity(ex('50=1X49='), 50) == pytest.approx(98.0)
    assert worst_window_identity(ex('50=1X49='), 100) == pytest.approx(99.0)
    assert worst_window_identity(ex('50=1X49='), 1000) == pytest.approx(99.0)
    assert worst_window_identity(ex('50=1X49='), 1000000) == pytest.approx(99.0)


def test_get_cigar():
    assert get_cigar('ACGATCGACATCACGACT', 'ACGATCGACATCACGACT', 'edlib') == '18='
    assert get_cigar('ACGATCGACATCACGACT', 'ACGATCGATATCACGACT', 'edlib') == '8=1X9='
