#!/usr/bin/env python3
"""Check for presence of a splint sequence in consensus reads using mappy.

Usage: check_splint_in_consensus.py -s SPLINT -r reads.fa [-n N] [-m MIN_MATCH_LEN]

Prints a summary of how many reads contain the splint, match lengths, positions,
and a few example alignments showing the read and splint with match markers.
"""

from __future__ import annotations

import argparse
import sys
import shutil
from collections import Counter
from typing import List, Tuple

import gzip
import re

# guarded optional imports: mappy is required at runtime; matplotlib is optional for plotting
try:
	import mappy as mp
except Exception:
	mp = None

try:
	import matplotlib.pyplot as plt
except Exception:
	plt = None

	# We use the list-style cigar attribute exposed by the installed mappy (list of [len, op_int]).


def parse_args():
	p = argparse.ArgumentParser(description="Check for splint sequence in consensus reads using mappy")
	p.add_argument("-s", "--splint", required=True, help="Splint sequence (DNA)")
	p.add_argument("-r", "--reads", required=True, help="FASTA/FASTA.GZ file of consensus reads")
	p.add_argument("-n", "--examples", type=int, default=5, help="Number of example alignments to print")
	p.add_argument("-m", "--min-match-len", type=int, default=8, help="Minimum match length to consider as splint present")
	p.add_argument("--min-id", type=float, default=0.0, help="Minimum percent identity (0-100) to consider a valid hit")
	p.add_argument("--no-rc", action="store_true", help="Do not search reverse-complement of splint")
	p.add_argument("--plot", help="Write coverage plot to this PNG filename (requires matplotlib)")
	p.add_argument("--scale", action='store_true', help="Scale read coordinates to 0-1 on x-axis")
	p.add_argument("--scale-bins", type=int, default=200, help="When --scale: number of bins across normalized read (default 200)")
	p.add_argument("--bin-size", type=int, default=0, help="Bin coverage by match length; 0 = no binning, otherwise bin size in bases")
	p.add_argument("--all-hits", action='store_true', help="Report all hits per read that meet thresholds (default: best only)")
	p.add_argument("--max-hits", type=int, default=5, help="When --all-hits is set, limit to this many hits per read")
	return p.parse_args()


def read_fasta_simple(path: str):
	"""A minimal FASTA reader that yields (header, seq). Works with plain FASTA only."""
	name = None
	seq_lines: List[str] = []
	open_fn = open
	if path.endswith('.gz'):
		open_fn = gzip.open
	with open_fn(path, 'rt') as fh:
		for line in fh:
			line = line.rstrip('\n')
			if not line:
				continue
			if line.startswith('>'):
				if name is not None:
					yield name, ''.join(seq_lines)
				name = line[1:].split()[0]
				seq_lines = []
			else:
				seq_lines.append(line)
		if name is not None:
			yield name, ''.join(seq_lines)


def rc(seq: str) -> str:
	comp = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
	return seq.translate(comp)[::-1]


# IUPAC ambiguity map: maps splint/base letter to set of matching bases
IUPAC = {
	'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'),
	'R': set('AG'), 'Y': set('CT'), 'S': set('GC'), 'W': set('AT'),
	'K': set('GT'), 'M': set('AC'), 'B': set('CGT'), 'D': set('AGT'),
	'H': set('ACT'), 'V': set('ACG'), 'N': set('ACGT')
}


def iupac_match(read_base: str, splint_base: str) -> bool:
	"""Return True if read_base matches splint_base considering IUPAC ambiguity on splint_base.
	splint_base may be lowercase or uppercase; read_base likewise.
	"""
	if not read_base or not splint_base:
		return False
	a = read_base.upper()
	b = splint_base.upper()
	# if splint base is standard ACGT, quick check
	if b in IUPAC:
		return a in IUPAC[b]
	# unknown character: fall back to direct equality
	return a == b


def format_alignment(read_seq: str, splint_seq: str, rstart: int, rend: int, qstart: int, qend: int) -> Tuple[str, str, str]:
	"""Return strings (read, splint_aligned, marker) where marker shows matches with '|' or '-' for mismatches.
	This function aligns the substring positions directly (no gaps) — good for rough visualization when the
	match reported by mappy is ungapped. If gapped alignment is desired, examples will still be useful.
	"""
	# extract parts
	r_sub = read_seq[rstart:rend]
	q_sub = splint_seq[qstart:qend]
	marker = []
	out_q = []
	out_r = []
	for a, b in zip(r_sub, q_sub):
		out_r.append(a)
		out_q.append(b)
		if iupac_match(a, b):
			marker.append('|')
		else:
			marker.append('-')
	return ''.join(out_r), ''.join(out_q), ''.join(marker)


def parse_cigar(cigar) -> List[Tuple[str, int]]:
	"""Normalize cigar information to a list of (op_char, length).
	Accepts either a mappy-style list of tuples (op_int, len) or a cigar string like '10M1I5D'.
	"""
	if cigar is None:
		return []
	# assume list-style cigar is list of [len, op_int]
	if isinstance(cigar, (list, tuple)) and len(cigar) and isinstance(cigar[0], (list, tuple)):
		op_map = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
		return [(op_map.get(op, 'M'), ln) for ln, op in cigar]
	# else assume string like '10M1I5D'
	parts = re.findall(r'(\d+)([MIDNSHP=X])', str(cigar))
	return [(op, int(length)) for length, op in parts]


def build_gapped_alignment(hit, read_seq: str, splint_seq: str, strand: str):
	"""Build gapped alignment strings from hit cigar. Returns (r_aln, s_aln, marker, qstart,qend,rstart,rend, matches, aln_len)
	If cigar not available, falls back to ungapped substring comparison.
	"""
	# get positions
	try:
		qpos = hit.q_st
		rpos = hit.r_st
	except Exception:
		# fallback: use 0
		qpos = 0
		rpos = 0

	# Use the mappy `cigar` attribute directly (installed mappy exposes a list-style cigar
	# and also provides `cigar_str`), keep it simple and rely on our parse_cigar() to handle
	# list or string formats.
	cigar = getattr(hit, 'cigar', None)

	ops = parse_cigar(cigar)
	if not ops:
		# fallback to ungapped
		try:
			qstart, qend = hit.q_st, hit.q_en
			rstart, rend = hit.r_st, hit.r_en
			r_aln, s_aln, marker = format_alignment(read_seq, splint_seq, rstart, rend, qstart, qend)
			matches = marker.count('|')
			return r_aln, s_aln, marker, qstart, qend, rstart, rend, matches, len(marker)
		except Exception:
			return '', '', '', 0, 0, 0, 0, 0, 0

	r_aln = []
	s_aln = []
	marker = []
	qpos0 = qpos
	rpos0 = rpos
	matches = 0
	aln_len = 0
	for op, ln in ops:
		if op in ('M', '=', 'X'):
			# consume both
			r_block = read_seq[qpos:qpos+ln]
			s_block = splint_seq[rpos:rpos+ln]
			for a, b in zip(r_block, s_block):
				r_aln.append(a)
				s_aln.append(b)
				if iupac_match(a, b):
					marker.append('|')
					matches += 1
				else:
					marker.append('-')
			qpos += ln
			rpos += ln
			aln_len += ln
		elif op == 'I' or op == 'S':
			# insertion to reference; present in read only
			r_block = read_seq[qpos:qpos+ln]
			r_aln.extend(list(r_block))
			s_aln.extend(['-'] * ln)
			marker.extend([' '] * ln)
			qpos += ln
			aln_len += ln
		elif op == 'D' or op == 'N':
			# deletion from read; present in ref only
			s_block = splint_seq[rpos:rpos+ln]
			r_aln.extend(['-'] * ln)
			s_aln.extend(list(s_block))
			marker.extend([' '] * ln)
			rpos += ln
			aln_len += ln
		else:
			# unknown op: skip
			pass

	qstart = qpos0
	qend = qpos
	rstart = rpos0
	rend = rpos
	return ''.join(r_aln), ''.join(s_aln), ''.join(marker), qstart, qend, rstart, rend, matches, aln_len


def create_aligners(splint: str, splint_rc: str, no_rc: bool):
	"""Create and return (a, a_rc) aligners for splint and its rc (a_rc may be None)."""
	if mp is None:
		raise RuntimeError('mappy (mp) module is not available')
	a = mp.Aligner(seq=splint, preset='map-ont')
	if not a:
		a = mp.Aligner()
		a.add_reference(splint, name='splint')

	a_rc = None
	if not no_rc:
		a_rc = mp.Aligner(seq=splint_rc, preset='map-ont')
		if not a_rc:
			a_rc = mp.Aligner()
			a_rc.add_reference(splint_rc, name='splint_rc')
	return a, a_rc


# (Removed detect_cigar_mode - we rely on list-style cigar [len, op_int])



def select_accepted_hits(candidates, seq: str, splint: str, splint_rc: str, args) -> List[tuple]:
	"""Given mapped candidates, build alignments, filter by thresholds and return accepted hits.

	Each returned tuple matches the shape used later in main():
	(hit, strand, aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk)
	"""
	passing = []
	best_score = (-1.0, 0)
	best_hit = None
	best_metrics = None

	for hit, strand, aligner in candidates:
		# attempt to build alignment, falling back to ungapped
		try:
			r_aln, s_aln, mrk, qstart, qend, rstart, rend, matches, aln_len = build_gapped_alignment(
				hit, seq, (splint if strand == '+' else splint_rc), strand
			)
		except Exception as e:
			# If cigar parsing raised an explicit ValueError (ambiguous/unsafe), skip this hit to avoid
			# producing incorrect alignments.
			if isinstance(e, ValueError):
				print(f"Skipping hit due to cigar parse error: {e}", file=sys.stderr)
				continue
			# otherwise, try a conservative ungapped fallback using reported coordinates
			try:
				qstart, qend = hit.q_st, hit.q_en
				rstart, rend = hit.r_st, hit.r_en
				# format_alignment expects (read_seq, splint_seq, rstart, rend, qstart, qend)
				r_aln, s_aln, mrk = format_alignment(seq, (splint if strand == '+' else splint_rc), rstart, rend, qstart, qend)
				matches = mrk.count('|')
				aln_len = len(mrk)
			except Exception:
				# give up on this candidate
				continue

		if aln_len == 0:
			pct_id = 0.0
		else:
			pct_id = matches / aln_len * 100.0

		# thresholds
		if aln_len < args.min_match_len:
			continue
		if pct_id < args.min_id:
			continue

		passing.append((pct_id, matches, aln_len, hit, strand, qstart, qend, rstart, rend, r_aln, s_aln, mrk))
		score = (pct_id, matches)
		if score > best_score:
			best_score = score
			best_hit = (hit, strand)
			best_metrics = (aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk)

	accepted = []
	if args.all_hits:
		passing.sort(key=lambda x: (x[0], x[1]), reverse=True)
		for tup in passing[: args.max_hits]:
			pct_id, matches, aln_len, hit, strand, qstart, qend, rstart, rend, r_aln, s_aln, mrk = tup
			accepted.append((hit, strand, aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk))
	else:
		if best_hit is not None:
			aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk = best_metrics
			accepted.append((best_hit[0], best_hit[1], aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk))

	return accepted


def print_summary_and_examples(total, with_splint, match_lengths, positions, examples, args):
	print(f"Total reads: {total}")
	print(f"Reads with splint (min_len={args.min_match_len}, min_id={args.min_id}%): {with_splint}")
	if match_lengths:
		c = Counter(match_lengths)
		most_common = c.most_common(5)
		print("Match length counts (top 5):")
		for ln, cnt in most_common:
			print(f"  {ln}: {cnt}")
	if positions:
		starts = [p[1] for p in positions]
		ends = [p[2] for p in positions]
		print(f"Match start in read: min={min(starts)}, max={max(starts)}, median~={sorted(starts)[len(starts)//2]}")
		print(f"Match end in read: min={min(ends)}, max={max(ends)}, median~={sorted(ends)[len(ends)//2]}")

	if examples:
		print('\nExamples:')
		term_width = shutil.get_terminal_size((80, 20)).columns
		for ex in examples:
			(name, r_sub, q_sub, marker, strand, qstart, qend, rstart, rend, matches, aln_len, pct_id) = ex
			read_label = f"R[{qstart}:{qend}]:"
			spl_label = f"S[{rstart}:{rend}]:"
			label_w = max(len(read_label), len(spl_label)) + 1
			read_label = read_label.ljust(label_w)
			spl_label = spl_label.ljust(label_w)

			avail = max(20, term_width - label_w - 1)
			seq_len = max(len(r_sub), len(q_sub), len(marker))
			r_full = r_sub.ljust(seq_len)
			s_full = q_sub.ljust(seq_len)
			m_full = marker.ljust(seq_len)

			print(f"Name: {name}")
			for i in range(0, seq_len, avail):
				r_chunk = r_full[i:i+avail]
				m_chunk = m_full[i:i+avail]
				s_chunk = s_full[i:i+avail]
				rl = read_label if i == 0 else ' ' * label_w
				sl = spl_label if i == 0 else ' ' * label_w
				print(f"{rl}{r_chunk}")
				print(f"{' ' * (label_w)}{m_chunk}")
				print(f"{sl}{s_chunk}")
			print(f"Identity: {matches}/{aln_len} ({pct_id:.1f}%)\n")


def plot_coverage(args, cov_read, cov_scaled, bins, bins_scaled, with_splint):
	if not args.plot:
		return
	if plt is None:
		print('matplotlib required for --plot but not importable', file=sys.stderr)
		return

	if args.scale:
		if cov_scaled is None:
			print('No scaled read-coordinate coverage to plot')
			return
		x = [i / float(len(cov_scaled) - 1) if len(cov_scaled) > 1 else 0.0 for i in range(len(cov_scaled))]
		plt.figure(figsize=(10, 3))
		if args.bin_size and bins_scaled:
			for idx, b in enumerate(bins_scaled):
				plt.plot(x, b[:len(x)], lw=1, label=f'bin{idx*args.bin_size+1}-{(idx+1)*args.bin_size}')
			plt.legend()
		else:
			plt.plot(x, cov_scaled, lw=1)
		plt.xlabel('read position (scaled 0-1)')
	else:
		if not cov_read:
			print('No read-coordinate coverage to plot')
			return
		x = list(range(len(cov_read)))
		plt.figure(figsize=(10, 3))
		if args.bin_size and bins:
			for idx, b in enumerate(bins):
				plt.plot(x, b[:len(x)], lw=1, label=f'bin{idx*args.bin_size+1}-{(idx+1)*args.bin_size}')
			plt.legend()
		else:
			plt.plot(x, cov_read, lw=1)
		plt.xlabel('read position (bases)')

	plt.ylabel('coverage (accepted hits)')
	plt.title(f'Read-coordinate coverage (accepted_hits={with_splint})')
	plt.tight_layout()
	plt.savefig(args.plot)
	print(f'Wrote coverage plot to {args.plot}')


def main():
	args = parse_args()

	# ensure mappy was importable at module import time
	if mp is None:
		print('Error: mappy is required but not importable', file=sys.stderr)
		sys.exit(2)

	splint = args.splint
	splint_rc = rc(splint)

	# create aligners for splint and its reverse complement
	a, a_rc = create_aligners(splint, splint_rc, args.no_rc)

	# we rely on the list-style cigar attribute exposed by mappy (list of [len, op_int])

	total = 0
	with_splint = 0
	match_lengths = []
	positions = []
	examples = []
	# per-base coverage across splint (use splint length)
	splint_len = len(splint)
	cov = [0] * splint_len
	# per-base coverage across reads (dynamic length)
	cov_read = []
	max_read_pos = 0
	# scaled coverage across normalized read coordinates (fixed bins) when args.scale
	cov_scaled = None
	bins_scaled = None
	# binned coverage by match length
	bins = None

	# iterate reads
	for name, seq in read_fasta_simple(args.reads):
		total += 1
		# collect candidate hits (forward and rc)
		candidates = []
		for h in a.map(seq):
			candidates.append((h, '+', a))
		if a_rc is not None:
			for h in a_rc.map(seq):
				candidates.append((h, '-', a_rc))

		# evaluate and select accepted hits using helper
		accepted = select_accepted_hits(candidates, seq, splint, splint_rc, args)
		# record accepted hits
		for hit_obj, strand, aln_len, matches, pct_id, qstart, qend, rstart, rend, r_aln, s_aln, mrk in accepted:
				with_splint += 1
				match_lengths.append(aln_len)
				positions.append((name, rstart, rend, qstart, qend, strand))
				# update per-base coverage on splint (rstart:rend)
				if rstart < rend and rstart >= 0:
					rs = max(0, rstart)
					re = min(splint_len, rend)
					for i in range(rs, re):
						cov[i] += 1
				# update per-base coverage on read (qstart:qend)
				if qstart is not None and qend is not None and qstart < qend:
					# if not scaling, accumulate absolute read-coordinate coverage as before
					if not args.scale:
						if qend > max_read_pos:
							# extend cov_read
							cov_read.extend([0] * (qend - max_read_pos))
							max_read_pos = qend
						for i in range(max(0, qstart), qend):
							cov_read[i] += 1
						# update bins if requested (absolute coordinates)
						if args.bin_size and args.bin_size > 0:
							bin_idx = max(0, (aln_len - 1) // args.bin_size)
							if bins is None:
								bins = []
							# create bins up to current index
							while len(bins) <= bin_idx:
								bins.append([0] * max_read_pos)
							# ensure existing bins are long enough
							for b in bins:
								if len(b) < max_read_pos:
									b.extend([0] * (max_read_pos - len(b)))
							# increment bin coverage
							for i in range(max(0, qstart), qend):
								bins[bin_idx][i] += 1
					else:
						# scaling: accumulate into normalized bins per read
						if cov_scaled is None:
							cov_scaled = [0] * args.scale_bins
							bins_scaled = [] if not args.bin_size else []
						read_len = len(seq)
						if read_len <= 1:
							# map everything to bin 0
							for i in range(max(0, qstart), qend):
								cov_scaled[0] += 1
						else:
							denom = float(read_len - 1)
							for i in range(max(0, qstart), qend):
								frac = i / denom
								bidx = int(frac * (args.scale_bins - 1))
								cov_scaled[bidx] += 1
							# update scaled bins by match length if requested
							if args.bin_size and args.bin_size > 0:
								bin_idx = max(0, (aln_len - 1) // args.bin_size)
								# ensure bins_scaled has enough rows
								while len(bins_scaled) <= bin_idx:
									bins_scaled.append([0] * args.scale_bins)
								for i in range(max(0, qstart), qend):
									frac = i / denom if read_len > 1 else 0.0
									bidx = int(frac * (args.scale_bins - 1))
									bins_scaled[bin_idx][bidx] += 1
				if len(examples) < args.examples:
					max_display = 200
					r_disp = r_aln
					s_disp = s_aln
					mrk_disp = mrk
					if len(r_disp) > max_display:
						r_disp = r_disp[:90] + '...' + r_disp[-90:]
						s_disp = s_disp[:90] + '...' + s_disp[-90:]
						mrk_disp = mrk_disp[:90] + '...' + mrk_disp[-90:]
					examples.append((name, r_disp, s_disp, mrk_disp, strand, qstart, qend, rstart, rend, matches, aln_len, pct_id))

	print_summary_and_examples(total, with_splint, match_lengths, positions, examples, args)

	plot_coverage(args, cov_read, cov_scaled, bins, bins_scaled, with_splint)


if __name__ == '__main__':
	main()
