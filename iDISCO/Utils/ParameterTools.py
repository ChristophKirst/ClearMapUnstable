# -*- coding: utf-8 -*-

"""
ParameterTools

Notes:
    - provides simple formatting tools to handle / print parameter organized as key:value paris


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

def writeParameter(head = None, out = None, **args):
    if head is None:
        head = '';
        
    keys = args.keys();
    vals = args.values();
    parsize = max([len(x) for x in keys]);
    
    s = [head + ' ' + keys[i].ljust(parsize) + ': ' + str(vals[i]) for i in range(len(keys))];
    
    if out is None:
        return '\n'.join(s)
    else:
        out.write('\n'.join(s));

    
def joinParameter(*args):
    """Joins dicts values are defined by first key : value pair if key has multiple occurences"""
    
    keyList = [x.keys() for x in args];
    n = len(args);
    
    keys = [];
    values = [];
    for i in range(n):
        values = values + [args[i][k] for k in keyList[i] if k not in keys];
        keys   = keys + [k for k in keyList[i] if k not in keys];
    
    return {keys[i] : values[i] for i in range(len(keys))}
    
    
    
    