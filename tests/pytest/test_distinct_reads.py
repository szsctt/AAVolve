import gzip
import os
import subprocess


def write_tsv(path, lines, gz=False):
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write('\n'.join(lines) + '\n')
    else:
        with open(path, 'wt') as f:
            f.write('\n'.join(lines) + '\n')


def read_tsv_lines(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as f:
            return [l.rstrip('\n') for l in f]
    else:
        with open(path, 'rt') as f:
            return [l.rstrip('\n') for l in f]


def test_single_column_assigned(tmp_path):
    # header only + several read ids (single-column case)
    inp = tmp_path / "assigned.tsv"
    out_counts = tmp_path / "counts.tsv"
    out_members = tmp_path / "members.tsv"

    lines = ["read_id", "r1", "r2", "r3"]
    write_tsv(str(inp), lines)

    cmd = ["python3", "aavolve/distinct_reads.py", "-i", str(inp), "-o", str(out_counts), "-m", str(out_members)]
    subprocess.check_call(cmd)

    counts_lines = read_tsv_lines(str(out_counts))
    members_lines = read_tsv_lines(str(out_members))

    assert counts_lines[0] == 'count'
    assert counts_lines[1].strip() == '3'

    # members should map G1 to r1,r2,r3
    assert members_lines[0] == 'group_id\tread_id'
    assert set(members_lines[1:]) == {'G1\tr1', 'G1\tr2', 'G1\tr3'}


def test_multi_column_assigned(tmp_path):
    # multi-column case: read_id + two parent columns
    inp = tmp_path / "assigned2.tsv.gz"
    out_counts = tmp_path / "counts2.tsv.gz"
    out_members = tmp_path / "members2.tsv.gz"

    # header: read_id, p1, p2
    lines = ['read_id\tp1\tp2', 'r1\ta\tb', 'r2\ta\tb', 'r3\tc\td']
    write_tsv(str(inp), lines, gz=True)

    cmd = ["python3", "aavolve/distinct_reads.py", "-i", str(inp), "-o", str(out_counts), "-m", str(out_members)]
    subprocess.check_call(cmd)

    counts_lines = read_tsv_lines(str(out_counts))
    members_lines = read_tsv_lines(str(out_members))

    # counts header should start with 'count'
    assert counts_lines[0].split('\t')[0] == 'count'

    # There should be two groups (one with count 2 and one with count 1)
    counts_vals = [int(l.split('\t')[0]) for l in counts_lines[1:]]
    assert sorted(counts_vals, reverse=True) == [2, 1]

    # members should include G1 and G2 mapping
    assert members_lines[0] == 'group_id\tread_id'
    groups = {l.split('\t')[0] for l in members_lines[1:]}
    assert groups == {'G1', 'G2'}
