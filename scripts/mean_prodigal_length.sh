#!/usr/bin/env bash

# This script runs Prodigal (github.com/hyattpd/Prodigal) to produce a mean protein length for an
# assembly. This allows for the relative ranking of alternative assemblies of the same genome.
# Higher mean protein lengths suggest a better assembly, as indel errors create frame-shifts which
# fragment coding sequences.

# This script takes one positional argument: assembly filename
# It produces one file: the assembly filename with ".prod" appended to the end

# Example usage:
# mean_prodigal_length.sh assembly.fasta

# Requirements: Prodigal and seqtk

# Copyright 2022 Ryan Wick (rrwick@gmail.com)

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You
# should have received a copy of the GNU General Public License along with this program. If not,
# see <https://www.gnu.org/licenses/>.

if [[ -z "$1" ]] ; then echo "you must provide an assembly filename"; exit 1; fi
if [[ ! -f "$1" ]] ; then echo $1" does not exist"; exit 1; fi

prod=$(mktemp)

prodigal -a "$prod" -i $1 > /dev/null
seqtk seq "$prod" | awk 'NR % 2 == 0' | awk '{sum += length($0); n++} END {print sum/n;}' > "$1".prod

rm "$prod"
