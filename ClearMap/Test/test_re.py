# -*- coding: utf-8 -*-
"""
Module to check and convert of regular and file expressions 

Utility module for source modules FileList and FileArray

"""
import sre_parse as sre

def isExpression(expression, groups = None, nPatterns = None, exclude = None, verbose = False):
  """Checks if the regular expression fullfill certain criteria
  
  Arguments:
    expression (str): regular exrpession to check
    groups (list or None): list of group names that should be present
    nPatterns (int or None): number of patterns to expect
    exclude (list or None): exculde these tokens
    verbose (bool): if True print reason for expression to not fullfil desired criteria
    
  Returns
    bool: True if the expression fullfills the desired criteria
  """
  
  #parse regular expression 
  p = sre.parse(source);

  #group patterns
  gd = p.pattern.groupdict
    
  if groups is None:
    groups = [];
  
  for gn in gd.keys():
    if gn not in groups:
      if verbose:
        print 'Expression does contain a non required group %s!' % gn;
      return False;
  
  for gn in groups:
    if gn not in gd.keys():
      if verbose:
        print 'Expression does not contain required group %s!' % gn;
      return False;

  if exclude is None:
    exclude = [];
  
  n = 0;
  for i,l in enumerate(p):
    lit = l[0];
    if lit != 'literal' and lit not in exclude:
      n = n + 1;
      if nPatterns is not None and n > nPatterns:
        if verbose:
          print 'Expression has more than %d regular expression patterns!' % nPatterns;
        return False;
    
  return True;

 
def test():
  


if __name__ == '__main__':
  test();