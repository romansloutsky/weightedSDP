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


#===============================================================================
# Substitution Matrix Class
#===============================================================================

class SubstitutionMatrix(object):
  '''
  Only identity comparison (1 if identical, 0 if not) is currently implemented.
  '''
  PROTEIN_LETTERS = frozenset(PROTEIN_ALPHABET)
  
  class KeyReturningDefaultDict(defaultdict):
    def __missing__(self,key):
      self[key] = self.default_factory(key)
      return self[key]
  
  def __init__(self,matrix=None):
    if matrix is None:
      self._matrix = self.KeyReturningDefaultDict(lambda x: len(x) % 2
                                                  if x < self.PROTEIN_LETTERS
                                                  else {}[x])
    else:
      raise NotImplementedError("Only identity comparison scoring is currently implemented")
  
  def __call__(self,item1,item2=None):
    '''
    The pair of residues to be compared may be passed as an iterable or as
    separate arguments
    '''
    score_this = frozenset(item1) if item2 is None else frozenset({item1,item2})
    try:
      return self._matrix[score_this]
    except KeyError:
      as_list = list(score_this - self.PROTEIN_LETTERS)
      if len(as_list) == 1:
        error_head = "Character "+as_list.pop()+" is "
      else:
        error_head = "Characters "+repr(as_list)+" are "
      raise KeyError(error_head+"not part of the protein alphabet "+\
                     repr(PROTEIN_ALPHABET))

  @classmethod
  def get_matrix(cls,matrix=None):
    '''
    Pass-through constructor only creates an instance when the argument is not
    already an instance
    '''
    if isinstance(matrix,cls):
      return matrix
    else:
      return cls(matrix)
