from collections.abc import Mapping, Callable
from copy import deepcopy
from io import StringIO
from math import pi
from numbers import Integral, Real
import os

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .reaction import Reaction
from .ace import Table, get_metadata, get_table
from .data import ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, get_head_record, get_tab1_record, get_list_record
from .function import Tabulated1D

_REACTION_NAME = {
  4: ('(alpha,n)', 'an'),
  16: ('(alpha,2n)', 'a2n'),
  17: ('(alpha,3n)', 'a3n'),
  37: ('(alpha,4n)', 'a4n'),
  152: ('(alpha,5n)', 'a5n'),
  153: ('(alpha,6n)', 'a6n'),
  152: ('(alpha,7n)', 'a7n'),
  153: ('(alpha,8n)', 'a8n'),
}

class IncidentAlpha(EqualityMixin):
  r"""Alpha particle interaction data

  This class stores data about alpha particle cross-section
  read from the ENDF-databases if possible
  """

  def __init__(self, atomic_number):
    self._atomic_number = atomic_number
    self.reactions = {}

  #REMOVEME these bits are copied in from photons.py
  def __contains__(self,mt):
    return mt in self.reactions

  def __getitem__(self, mt):
    if mt in self.reactions:
      return self.reactions[mt]
    else:
      raise KeyError('No reaction with MT={}.'.format(mt))

  def __iter__(self):
    return iter(self.reaction.values())

  @property
  def atomic_number(self):
    return self.atomic_number

  @atomic_number.setter
  def atomic_number(self, atomic_number):
    cv.check_type('atomic number', atomic_number, Integral)
    cv.check_greater_than('atomic number', atomic_number, 0, True)
    self._atomic_number = atomic_number

  @property
  def name(self):
    return ATOMIC_SYMBOL[self.atomic_number]

  @classmethod
  def from_endf(cls, ev_or_filename):
    """Generate incident alpha-particle data from ENDF evaluation
    Parameters
    ----------
    ev_or_filename : openmc.data.evaluation or filename
      of ENDF-file that contains the data
    """
    if isinstance(alpha, Evaluation):
      ev = ev_or_filename
    else:
      ev = Evaluation(ev_or_filename)

    Z = ev.target['atomic_number']
    data = cls(Z)

    # Read each reaction
    for mf, mt, nc, mod in ev.reaction_list:
      if mf == 3:
        data.reactions[mt] = Reaction.from_endf(ev, mt)

    data._evaluation = ev
    return data
