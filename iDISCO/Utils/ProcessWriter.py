# -*- coding: utf-8 -*-

"""
PrcessWriter

provides simple formatting tools to print text with process header


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys

class ProcessWriter(object):
    def __init__(self, process = 0):
        self.process = process;
    
    def writeString(self, text):
        return ("Process %d: " % self.process) + str(text);
    
    def write(self, text):
        print self.writeString(text);
        sys.stdout.flush();

    