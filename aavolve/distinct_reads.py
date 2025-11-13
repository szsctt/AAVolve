#!/usr/bin/env python3
"""
Create collapsed parent-counts and a members mapping from assigned parents TSV.

This replaces the previous bash implementation in the Snakemake rule with a
single, testable Python module that correctly handles gzipped inputs/outputs.

Behavior:
- If the assigned TSV has only a single column (read_id) the module writes:
  - counts file with header 'count' and a single row containing the total
    number of read IDs (aggregate)
  - members file with header 'group_id\tread_id' and each read mapped to
    a single group 'G1'
- Otherwise the module collapses identical rows (excluding the first column)
  and writes the counts file with header (first column replaced by 'count')
  and the members mapping that assigns synthetic group IDs (G1, G2, ...) to
  each unique group and lists the read IDs for that group.

The module is gz-aware and uses aavolve.utils.use_open for I/O.
"""

import argparse
import csv
from collections import defaultdict

from aavolve.utils import use_open
def main():
    parser = argparse.ArgumentParser(description='Collapse assigned parents to counts and members mapping')
    parser.add_argument('-i', '--input', required=True, help='Assigned parents TSV (may be gz)')
    parser.add_argument('-o', '--output', required=True, help='Output counts TSV (may be gz)')
    parser.add_argument('-m', '--members', required=True, help='Output members mapping TSV (may be gz)')
    args = parser.parse_args()

    infile = args.input
    out_counts = args.output
    out_members = args.members

    # Read header and determine number of columns
    with use_open(infile, 'rt') as fh:
        header_line = fh.readline()
        if header_line == '':
            raise ValueError(f"Input file '{infile}' is empty")
        header_line = header_line.rstrip('\n')
        headers = header_line.split('\t')

        # If only read_id present -> single-column case
        if len(headers) <= 1:
            # Count total non-empty read ids and write aggregate counts
            read_ids = []
            for line in fh:
                line = line.strip()
                if line:
                    read_ids.append(line)

            # write counts (single aggregated count)
            with use_open(out_counts, 'wt') as cf:
                cf.write('count\n')
                cf.write(str(len(read_ids)) + '\n')

            # write members mapping: group G1 -> all read ids
            with use_open(out_members, 'wt') as mf:
                mf.write('group_id\tread_id\n')
                for rid in read_ids:
                    mf.write(f'G1\t{rid}\n')

            return

        # Normal case: collapse by the remaining columns
        # We'll map from group_key -> list of read_ids
        group_map = defaultdict(list)
        # header for counts: replace first column name with 'count'
        counts_header = ['count'] + headers[1:]

        # read remaining lines
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row:
                continue  # skip empty lines
            # error if row is short
            assert len(row) == len(headers), f"Row has {len(row)} columns but expected {len(headers)}: {row}"
            read_id = row[0]
            group_key = '\t'.join(row[1:])
            assert group_key != "", f"Empty group key for read_id {read_id} in row: {row}"
            group_map[group_key].append(read_id)

    # create counts text sorted by descending count
    groups = sorted(group_map.items(), key=lambda kv: len(kv[1]), reverse=True)

    # write counts file
    with use_open(out_counts, 'wt') as cf:
        cf.write('\t'.join(counts_header) + '\n')
        for group_idx, (group_key, members) in enumerate(groups, start=1):
            count = len(members)
            cf.write(f"{count}\t{group_key}\n")

    # write members mapping with synthetic group ids
    with use_open(out_members, 'wt') as mf:
        mf.write('group_id\tread_id\n')
        for group_idx, (group_key, members) in enumerate(groups, start=1):
            gid = f'G{group_idx}'
            for rid in members:
                mf.write(f"{gid}\t{rid}\n")


if __name__ == '__main__':
    main()
