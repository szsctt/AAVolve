# this script extracts the variants in each read relative to a reference

#	Write to tab-separated output file
# 	Columns are:
# 	 - reference name
# 	 - position
# 	 - query name
# 	 - variant string (e.g. S492A, where S is reference and A is query)
# 	 - ref bases (. if insertion)
# 	 - variant bases (. if deletion)
#    - if variant changes an amino acid

# Alternative for smaller output file is with columns:
# 	 - query name
# 	 - position
# 	 - ref bases (. if insertion)
# 	 - variant bases (. if deletion)
#    - if variant changes an amino acid

# positions - always relative to reference:
# - substiutions: one-based position of substution
# - deletions: one-based positions of deleted bases (e.g. 50_60)
# - insertions: one-based positions of bases either side of inserted bases (e.g. 50_51)

# mutation strings:
# mutation numbering is one-based:
# - substituion of base 250 C (in reference) to A (in query): C250A
# - insertion of GC between bases 255 and 256: 255_256insGC
# - deletion of a single base, G that was at position 290: G290del
# - deletion of four bases, GTCC, which were positions 250-253: G250_C253del

import argparse
import os
from sys import argv

import pysam
import tqdm

from Bio.Seq import Seq
from Bio import SeqIO

from scripts.utils import Substitution, Deletion, Insertion
from scripts.utils import use_open

def main():

  args = parse_args(argv[1:])

  get_all_variants(args.i, args.r, args.o, args.must_start_before, args.must_end_after, args.smaller_output, args.aa_change_in_isolation)
  

def get_all_variants(samfile, reffile, outfile, must_start_before, must_end_after, smaller_output, aa_isolation):

  # check that index exists for samfile
  check_index_exists(samfile)

  ref_names = get_reference_names(samfile)

  ref_seqs = SeqIO.to_dict(SeqIO.parse(reffile, "fasta"))

  with use_open(outfile, "wt") as out:
    
    # write header
    write_header(out, smaller_output)

    # get variants for each read
    for read_id, var in get_variants(samfile, ref_seqs, must_start_before,
                                   must_end_after, aa_isolation):
        
        # write variants to file
        write_variant(read_id, var, ref_names, out, smaller_output)

def parse_args(argv):

  parser = argparse.ArgumentParser("Extract variants relative to reference for each read from sam file")
  parser.add_argument("-i", help="Sam file for input", required=True)
  parser.add_argument("-r", help="Reference file for input", required=True)
  parser.add_argument("-o", help="Tsv file for writing variants", default="variants.tsv")
  parser.add_argument('-b', "--must-start-before", type=int, default = 0, 
                      help="Reads will be excluded if they don't start "
                           "before this position in the reference")
  parser.add_argument('-a', "--must-end-after", type=int, default = -1, 
                      help="Reads will be excluded if they don't start "
                           "before this position in the reference "
                           "(-1 means end of reference)")
  parser.add_argument('-S', "--smaller-output", action="store_true", help="Output fewer columns")
  parser.add_argument('-I', "--aa-change-in-isolation", action="store_true", help="If set, substitutions are evaluated for amino acid changes in isolation, rather than in the context of the codon")
  args = parser.parse_args(argv)

  return args
  
def get_reference_names(samfile):
  """
  Get names of references for each read
  """
  
  refs = {}
  for read in pysam.AlignmentFile(samfile):
  
    refs[read.query_name] = read.reference_name
    
  return refs

def check_index_exists(samfile):
    """
    Check that index exists for samfile
    """
    if not os.path.exists(samfile + ".bai"):
        raise OSError(f"Index file {samfile}.bai not found")

def count_total_reads(pysam_handle):
    """
    Get the total number of mapped reads in the file 
    regardless of alignment flag (e.g. secondary alignments are included)
    """
    # count total reads for each reference sequence
    idx = pysam_handle.get_index_statistics()
    total = 0
    for i in idx:
        total += i.mapped

    return total

def get_variants(samfile, ref_seqs, start, end, aa_isolation):

  # open samfile
  reads = pysam.AlignmentFile(samfile)
  
  # get length of references
  if end == -1:
    refs = reads.lengths
  else:
    check_end = end
    
  discarded = 0
  started_aln = False
  
  # iterate over reads in samfile
  total = count_total_reads(reads) 
  for read in tqdm.tqdm(reads, desc="Processing reads", total=total):

    # skip unmapped reads
    if not read.is_mapped:
      continue
  
    # skip supplementary alignments
    if read.is_supplementary:
      continue
      
    # skip secondary alignments
    if read.is_secondary:
      continue
  
    # check that read covers at least from 'start' position in reference onwards:
    if start < read.reference_start:
      discarded += 1
      continue
    
    # check that read covers at least to 'end' position in reference
    if end == -1:
      check_end = refs[read.reference_id]  
    
    if read.reference_end < check_end:
      discarded += 1
      continue
    
    # list of mutation objects
    read_vars = []
    
    last_rpos, last_qpos = None, None
    started_aln = False

    offset = read.query_alignment_start - read.reference_start

    for qpos, rpos, rseq in read.get_aligned_pairs(with_seq=True):

    
      # skip deletions at the start of the read
      if not started_aln:
        if rpos is None:
          continue
        else:
          started_aln = True
      
      # keep track of position, for when there are deletions
      # doesn't get updated in insertions
      if rpos is not None:
        curr_refpos = rpos
      
      # substitution
      if qpos is not None and rpos is not None:
      
         # substituions are lower case
        if rseq.islower():

          # get info about variant
          query_bases = get_query_base(read, qpos)
          assert query_bases.upper() != rseq.upper()
          
          changes_aa = identify_aa_change(read, ref_seqs, qpos, rpos, offset, aa_isolation)

          sub = Substitution(rpos = rpos, rseq = rseq, 
                              qseq = query_bases, changes_aa=changes_aa)
        
          read_vars.append(sub)
        
      # deletion
      if qpos is None:
      
        # if no variants, add a new one
        if len(read_vars) == 0:
          read_vars.append(Deletion(rpos, last_qpos, rseq))
        
        # if the last mutation was not a deletion, add a new one
        elif read_vars[-1].var_type != "del":
          read_vars.append(Deletion(rpos, last_qpos, rseq))
          
        # if the last mutation was not the same deletion, add a new one
        elif read_vars[-1].last_qpos != last_qpos:
          read_vars.append(Deletion(rpos, last_qpos, rseq))
          
        # otherwise, update current deletion
        else:
          read_vars[-1].add_another_base(rseq)

        offset -= 1
        
      # insertion
      if rpos is None:
      
        # get query sequence at this position
        qbase = get_query_base(read, qpos)
        
        # if no other variants, add new insertion
        if len(read_vars) == 0:
          read_vars.append(Insertion(last_rpos, qpos, qbase))
          
        # if last mutation wasn't an insertion, add new insertion   
        elif read_vars[-1].var_type != "ins":
          read_vars.append(Insertion(last_rpos, qpos, qbase))
        
        # if last insertion wasn't the same as this one, add new insertion
        elif read_vars[-1].last_rpos != last_rpos:
          read_vars.append(Insertion(last_rpos, qpos, qbase))
        
        # otherwise, append to insertion
        else:
          read_vars[-1].add_another_base(qbase)
        
        offset += 1
      
      # update rpos if not in insertion
      if rpos is not None:
        last_rpos = rpos
       
      # update qpos if not in deletion 
      if qpos is not None:
        last_qpos = qpos
    
    yield read.query_name, read_vars
  print(f"Discarded {discarded} reads")

def identify_aa_change(read, ref_seqs, qpos, rpos, offset, aa_isolation):
    """
    Identify variants that result in changes to amino acids

    If aa_isolation is false, then variants are evaluated in the context of the codon 
    That is, if there are multiple subsitutions in the same codon, then they are evaluated together
    Otherwise variants are evaluated in isolation

    Offset is the caused by insertions and deletions from the query sequence
    that cause the read to be shifted relative to the reference
    """


    # get reference codon - asume first codon starts at posiiton 0 in reference
    rcodon_start = rpos // 3 * 3
    rcodon_end = rcodon_start + 3
    rcodon = str(ref_seqs[read.reference_name][rcodon_start:rcodon_end].seq)
    assert len(rcodon) == 3

    if aa_isolation:
      # subtitute query base for reference base
      pos = rpos  % 3
      qcodon = rcodon[:pos] + read.query_sequence[qpos] + rcodon[pos + 1:]

    else:

      # get query codon
      qcodon_start = (qpos - offset) // 3 * 3 + offset
      qcodon_end = qcodon_start + 3
      qcodon = read.query_sequence[qcodon_start:qcodon_end]

      assert len(qcodon) == 3
    
    # compare amino acids
    return Seq(rcodon).translate() != Seq(qcodon).translate()

def get_query_base(read, qpos):

  return read.query_sequence[qpos]
 
def write_header(filehandle, smaller = False):

    if smaller:
       header = ["query_name", "pos", "ref_bases", "query_bases", "aa_change"]
    else:
        header = ["reference_name", "pos", "query_name", "var", "ref_bases", "query_bases", "aa_change"]

    filehandle.write("\t".join(header) + "\n")

def write_variant(read_id, variants, ref_names, filehandle, smaller = False):

    for var in variants:

        if smaller:
            filehandle.write(var.print_line(read_id))
        else:
            filehandle.write(var.print_line(read_id, ref_names[read_id]))

if __name__ == "__main__":
  main()
