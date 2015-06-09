# -*- coding: utf-8 -*-

"""
Timer

provides simple timing tools


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import time


class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose;
        self.start();

    def start(self):
        self.time = time.time();
    
    def reset(self):
        self.time = time.time();
    
    def elapsedTime(self, head = None, asstring = True):
        t = time.time();
        
        if asstring:
            t = self.formatElapsedTime(t - self.time);
            if head != None:
                return head + ": elapsed time: " + t;
            else:
                return "Elapsed time: " + t;
        else:
            return t - self.time;
    
    def printElapsedTime(self, head = None):
        print self.elapsedTime(head = head);
    
    def formatElapsedTime(self, t):
        m, s = divmod(t, 60);
        h, m = divmod(m, 60);
        
        return "%d:%02d:%02d" % (h, m, s);

