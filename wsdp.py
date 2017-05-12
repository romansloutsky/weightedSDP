from __future__ import division
import random
import math
import itertools
from collections import Counter,defaultdict
try:
  from Bio import SeqIO
  from Bio.Data.IUPACData import protein_letters as PROTEIN_ALPHABET
except ImportError:
  raise ImportError("Failed to import necessary Biopython components. "\
                    "Please install Biopython.")
