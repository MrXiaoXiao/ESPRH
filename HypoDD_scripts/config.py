""" Configure file for hypoDD interface
"""
import os
import numpy as np

class Config(object):
  def __init__(self):

    # 1. format input
    self.fsta_in = 'input/HYPO.sta'
    self.fsta_out = 'input/station.dat'
    self.fpha_in = 'input/merge.pha'
    self.fpha_out = 'input/phase.dat'
    self.dep_corr = 5 # avoid air quake

    # 2. format output
    self.out_ctlg = 'output/indonesia.ctlg'
    self.out_pha = 'output/indonesia.pha'
    self.out_pha_all = 'output/indonesia_all.pha'

