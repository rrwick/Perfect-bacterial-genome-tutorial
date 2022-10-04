#!/usr/bin/env bash

# This script runs ALE (github.com/sc932/ALE) to produce an ALE score for an assembly. ALE scores
# allow for the relative ranking of alternative assemblies of the same genome. Higher ALE scores
# (i.e. smaller magnitude, closer to zero) suggest a better assembly.

# This script takes four positional arguments:
# 1) assembly filename
# 2) Illumina read filename (first in pair)
# 3) Illumina read filename (second in pair)
# 4) number of threads to use for alignment

# It produces one file: the assembly filename with ".ale" appended to the end

# Example usage:
# ale_score.sh assembly.fasta reads/illumina_1.fastq.gz reads/illumina_2.fastq.gz 16

# Requirements: BWA and ALE

# Copyright 2022 Ryan Wick (rrwick@gmail.com)

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You
# should have received a copy of the GNU General Public License along with this program. If not,
# see <https://www.gnu.org/licenses/>.

if [[ -z "$1" ]] ; then echo "you must provide an assembly filename"; exit 1; fi
if [[ -z "$2" ]] ; then echo "you must provide read filenames"; exit 1; fi
if [[ -z "$3" ]] ; then echo "you must provide read filenames"; exit 1; fi
if [[ -z "$4" ]] ; then echo "you must provide a thread count"; exit 1; fi

if [[ ! -f "$1" ]] ; then echo $1" does not exist"; exit 1; fi
if [[ ! -f "$2" ]] ; then echo $2" does not exist"; exit 1; fi
if [[ ! -f "$3" ]] ; then echo $3" does not exist"; exit 1; fi

index=$(mktemp)
sam_file=$(mktemp)
ale_file=$(mktemp)

bwa index -p "$index" "$1"
bwa mem -t "$4" "$index" "$2" "$3" > "$sam_file"
ALE "$sam_file" "$1" "$ale_file"
grep "# ALE_score: " "$ale_file" | sed 's/# ALE_score: //' > "$1".ale

rm "$index" "$index".* "$sam_file" "$ale_file"
