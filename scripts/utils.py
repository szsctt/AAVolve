import gzip
import csv
import os

MAX_SEQS = 1000

def use_open(filename, *args, **kwargs):
    """
    Use gzip.open if file is gzipped
    """
    
    # file extension is .gz, assume should be opened with gzip
    if filename.endswith('.gz'):
      return gzip.open(filename, *args, **kwargs)
    # if file exists, check first two bytes to see if it's gzipped
    elif os.path.isfile(filename):
      with open(filename, 'rb') as f:
        magic = f.read(2)
      if magic == b'\x1f\x8b':
        return gzip.open(filename, *args, **kwargs)
    
    # assume non-gzipped file
    return open(filename, *args, **kwargs)
    
def get_repeats_from_r2c2_name(name):

    try:
        rep = int(name.split('_')[-2])
    except:
        raise ValueError(f"Can't get number of repeats from read '{name}'")
    return rep

def seq_generator(handle):
    """Generator that yields a sequence and its name from a FASTA file"""
    name, seq = '', ''
    for line in handle:
        if line.startswith('>'):
            if name != '':
                yield name, seq
            name = line.strip()
            seq = ''
        else:
            seq = seq + line.strip()
    if name != '':
        yield name, seq

def read_variant_file(filename):

  with use_open(filename, 'rt', newline='') as f:
    reader = csv.DictReader(f, delimiter="\t")
    for line in reader:
      yield line

def get_variant_type(ref, alt):
  """
  Get type of variant based on row of file
  """
  # insertion relative to reference
  if ref == '.':
    return 'ins'
  
  # deletion from reference
  elif alt == '.': 
    return 'del'
  
  return 'sub'

def get_variant(row):
  """
  Return a variant (Insertion, Deletion, Substitution) from a row from a file
  """
  var_type = get_variant_type(row['ref_bases'], row['query_bases'])
  

  # get info for substitution
  if var_type == 'sub':
     
    rpos = int(row['pos']) 
    rseq = row['ref_bases']
    qseq = row['query_bases']
    changes_aa = row['aa_change'] == 'True'


    assert len(rseq) == 1
    assert len(qseq) == 1
    assert changes_aa is True or changes_aa is False
    assert rpos >=0

    var = Substitution(rpos, rseq, qseq, changes_aa)
    
  # get info for insertion 
  elif var_type == 'ins':
     
    last_rpos = int(row['pos']) - 1
    qpos = None # we don't have this information, so set to None
    bases = row['query_bases']
    var = Insertion(last_rpos, qpos, bases)
     
  # get info for deletion
  elif var_type == 'del':
     
    
    pos = row['pos'].split('_')
    start_rpos = int(pos[0])
    end_rpos = int(pos[1])
    bases = row['ref_bases']
    
    var = Deletion(start_rpos, None, bases)
    var.end_rpos = end_rpos

  # check variant strings are equal
  if 'var' in row:
    assert str(var) == row['var']
   
  return var

def get_variants_set(filename):
  """
  Get a set of variants from file
  """
  vars = set()
  # read as tsv
  print(f"Reading variants from {filename}")
  
  for row in read_variant_file(filename):

    # get id for this variant
    var = get_variant(row)
    vars.add(var)

  return vars

def get_header(filename):
  """
  Get header from file
  """
  with use_open(filename, 'rt') as file_handle:
    header = file_handle.readline()

  if header == '':
    raise ValueError(f"File '{filename}' is empty")
  return header

def get_reference_name(filename):
  """
  Get reference name from file
  """

  # check if file is empty (no header or variants)
  with use_open(filename, 'rt') as file_handle:
    first_line = file_handle.readline()
    if first_line == '':
      raise ValueError(f"File '{filename}' is empty")

  # read first line to get reference name
  try:
    for row in read_variant_file(filename):
      ref = row['reference_name']
      return ref
  # no variants in file
  except StopIteration:
    pass
  # reference isn't output for shorter format 
  except KeyError:
    pass
  
  return None

def sort_var_names(variants):

    return sorted(variants, key=lambda x: int(x.split(':')[0].split('_')[0]))

def get_parents(parents_file):
    """
    Read in parents file
    Collect variants of the same type / position using variant id
    which consists of {pos}:{type}

    Return a dict of { id: {parent_name:var}}
    """
    parents = {}
    for row in read_variant_file(parents_file):
        # get variant from row  
        var = get_variant(row)
        var_id = var.var_id()

        # check if we've seen this variant before
        if var_id not in parents:
            parents[var_id] = {}

        # add variant to dict
        parent_name = row['query_name']
        parents[var_id][parent_name] = var

    return parents

def count_lines(filename):
  """
  Count number of lines in a file
  """
  count = 0
  with use_open(filename, 'rt') as file_handle:
    for line in file_handle:
      if line.strip() != '':
        count += 1
  return count
 
# internal representation of a subsitution
class Substitution:
  """
  Internal representation of a substitution
  """

  def __init__(self, rpos, rseq, qseq, changes_aa):
    """
    Initialise a substitution
    """
  
    self.var_type = "sub"
  
    # position in reference - 0-based numbering is used internally
    self.rpos = rpos
    
    # reference base
    self.rseq = rseq.upper()
    
    # query base
    self.qseq = qseq.upper()
    
    self.changes_aa = changes_aa

  def zero_pos(self):
    """
    Zero-based position of substitution
    """
    
    return self.rpos
    
  def one_pos(self):
    """
    One-based position of substitution
    """
    
    return self.rpos + 1

  def refbases(self):
    """
    Bases from the reference
    """
    
    return self.rseq
    
  def qbases(self):
    """
    Bases from the query
    """
    
    return self.qseq
  
  def print_line(self, query_name=None, ref_name=None):
    """
    Line for printing to file / stdout
    """
    assert query_name is not None

    if ref_name is None :
      fields = [query_name, self.zero_pos(), self.refbases(), self.qbases(), self.changes_aa]
      
    else:
      fields = [ref_name, self.zero_pos(), query_name, str(self), self.refbases(), self.qbases(), self.changes_aa]

    fields = [str(i) for i in fields]
    return "\t".join(fields) + "\n"
  
  def header(self, shorter=True):
    """
    Header for printing to file / stdout
    """

    if shorter:
      fields = ["query_name", "pos", "ref_bases", "query_bases", "aa_change"]
    
    else:
      fields = ["reference_name", "pos", "query_name", "var", "ref_bases", "query_bases", "aa_change"]
      
    return "\t".join(fields) + "\n"
    
  def var_id(self):
    """
    ID for variant consisting of position and type
    """

    return f"{self.zero_pos()}:{self.var_type}"

  def __str__(self):
    """
    Terse, one-based variant string
    """
    
    return f"{self.rseq}{self.rpos+1}{self.qseq}"


  def __repr__(self):
    """
    Verbose, zero-based information about variant
    """ 
  
    return (f"Substitution of query base '{self.qseq}' for reference base '{self.rseq}' "
            f"at 0-based reference position {self.rpos}")
  
  def __eq__(self, other):
    """
    Test for equality
    """

    return str(self) == str(other)
  
  def __hash__(self):
    """
    Hash of string
    """
    
    return hash(str(self))
            
# internal representation of an insertion     
class Insertion(Substitution):
  """
  Insertion is *into the reference* - i.e. query has bases that aren't in reference

  Example of an  insertion of more than one base (one-based var string: 50_51insCG)
     rpos  qpos  rseq  qseq
     49    60    A     A      <- last_rpos comes from this row
     None  61    None  C      <- start_qpos comes from this row
     None  62    None  G      
     50    63    T     T      <- end_qpos comes from this row
  """
  

  def __init__(self, last_rpos, qpos, first_query_base):
    """
    Initialise first base of insertion
    """
  
    self.var_type = "ins"
  
    # zero-based position of last reference base before insertion
    self.last_rpos = last_rpos
    
    # zero-based position of first query base in insertion
    self.start_qpos = qpos
    
    # zero-based position of first query base after insertion
    if qpos is not None:
      self.end_qpos = qpos + 1
    else:
      self.end_qpos = None
    
    # query bases that were inserted
    self.bases = first_query_base
    
    # insertions always change the amino acid sequence
    self.changes_aa = True
    
   
  def add_another_base(self, base):
    """
    Add another base to insertion
    """
  
    self.bases = self.bases + base
    self.end_qpos += 1
   
  def zero_pos(self):
    """
    Zero-based position of insertion - position after last reference base
    """
    
    return self.last_rpos + 1

  def one_pos(self):
    """
    One-based position of reference bases either side of insertion
    """
    
    pos = (self.last_rpos + 1, self.last_rpos + 2)
    return "_".join([str(i) for i in pos])

  def refbases(self):
    """
    Bases from the reference
    """
    
    return "."
    
  def qbases(self):
    """
    Bases from the query
    """
    
    return self.bases

    
  def __str__(self):
    """
    Terse, one-based variant string
    """
  
    return f"{self.last_rpos + 1}_{self.last_rpos + 2}ins{self.bases}"
    
  def __repr__(self):
    """
    Verbose, zero-based information about variant
    """
  
    return (f"Insertion of query bases {self.start_qpos}:{self.end_qpos} "
            f" ({self.bases}) after reference position {self.last_rpos}")
    
class Deletion(Substitution):
  """
  Deletion is *from the reference* - i.e. reference has bases that aren't in query

  Example of a deletion of one base (one-based var string: C51del)
     rpos  qpos  rseq  qseq
     49    60    A     A      <- last_qpos comes from this row
     50    None  C     None   <- start_rpos comes from this row
     51    61    T     T      <- end_rpos comes from this row

  Example of a deletion of more than one base (one-based var string: C51_G52del)
     rpos  qpos  rseq  qseq
     49    60    A     A      <- last_qpos comes from this row
     50    None  C     None   <- start_rpos comes from this row
     51    None  G     None   
     52    61    T     T      <- end_rpos comes from this row

  """

  def __init__(self, rpos, last_qpos, first_ref_base):
    """
    Initialise first base of deletion 
    """
  
    self.var_type = "del"
    
    # zero-based reference position of first base in deletion
    self.start_rpos = rpos
    
    # zero-based reference position of first base after deletion
    self.end_rpos = rpos + 1
    
    # zero-based query position of last base before deletion
    self.last_qpos = last_qpos
    
    # reference bases that were deleted
    self.bases = first_ref_base
    
    # deletions always change the amino acid sequence
    self.changes_aa = True
    
  def add_another_base(self, base):
    """
    Add another base to deletion
    """  
    self.bases = self.bases + base
    self.end_rpos += 1

  def zero_pos(self):
    """
    Zero-based positions of deleted bases
    """
    
    return f"{self.start_rpos}_{self.end_rpos}"

  def one_pos(self):
    """
    One-based positions of deleted bases
    """
    
    pos = (self.start_rpos + 1, self.end_rpos)
    return "_".join([str(i) for i in pos])

  def refbases(self):
    """
    Bases from the reference
    """
    
    return self.bases
    
  def qbases(self):
    """
    Bases from the query
    """
    
    return "."
    
    
  def __str__(self):
    """
    Terse, one-based variant string
     - deletion of a single base, G that was at position 290: G290del
     - deletion of four bases, GTCC, which were positions 250-253: G250_C253del  
    """
    
    if len(self.bases) == 1:
      return f"{self.bases}{self.start_rpos+1}del"
      
    return f"{self.bases[0]}{self.start_rpos+1}_{self.bases[-1]}{self.end_rpos}del"
    
  def __repr__(self):
    """
    Verbose, zero-based information about variant
    """
  
    return (f"Deletion of reference bases {self.start_rpos}:{self.end_rpos} "
            f" ({self.bases})")    
            
