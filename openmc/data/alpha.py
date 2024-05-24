from collections.abc import Mapping, Callable
from copy import deepcopy
from io import StringIO
from math import pi, sqrt
from numbers import Integral, Real
import os

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc import Material
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .reaction import Reaction
from .ace import Table, get_metadata, get_table
from .data import ATOMIC_NUMBER, ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, get_head_record, get_tab1_record, get_list_record
from .function import Tabulated1D

_REACTION_NAME = {
<<<<<<< HEAD
    4: ("(alpha,n)", "an"),
    16: ("(alpha,2n)", "a2n"),
    17: ("(alpha,3n)", "a3n"),
    37: ("(alpha,4n)", "a4n"),
    152: ("(alpha,5n)", "a5n"),
    153: ("(alpha,6n)", "a6n"),
    160: ("(alpha,7n)", "a7n"),
    161: ("(alpha,8n)", "a8n"),
}


def _init_alpha_cross_sections():
    an_cross_sections = {}
    # read the xml-file pointed to by config/environment var
    # to find the cross section datafiles
    import lxml.etree as ET

    try:
        fname = os.environ["OPENMC_AN_CROSS_SECTIONS"]
    except KeyError:
        print(
            'Error: Must specify a cross-section file in the environment var "OPENMC_AN_CROSS_SECTIONS" to use alpha-particle cross sections'
        )
        raise

    tree = ET.parse(fname)
    for xs in tree.iter():
        if xs.tag == "library" and "type" in xs.keys():
            if xs.attrib["type"] == "alpha_n":
                atom_symbol = xs.attrib["materials"]
                path = os.path.join(os.path.dirname(fname), xs.attrib["path"])
                an_cross_sections[atom_symbol] = path
    return an_cross_sections


_OPENMC_AN_CROSS_SECTIONS = _init_alpha_cross_sections()


class IncidentAlpha(EqualityMixin):
    r"""Alpha particle interaction data

    This class stores data about alpha particle cross-section
    read from the ENDF-databases if possible
    """

    def __init__(self, atomic_number):
        self._atomic_number = atomic_number
        self.reactions = {}

    # REMOVEME these bits are copied in from photons.py
    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            raise KeyError("No reaction with MT={}.".format(mt))

    def __iter__(self):
        return iter(self.reaction.values())

    @property
    def atomic_number(self):
        return self.atomic_number

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type("atomic number", atomic_number, Integral)
        cv.check_greater_than("atomic number", atomic_number, 0, True)
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
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        Z = ev.target["atomic_number"]
        data = cls(Z)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 3:
                data.reactions[mt] = Reaction.from_endf(ev, mt)

        data._evaluation = ev
        return data


def alpha_stopping_power(material: Material = None, E: float = 0.0) -> float:
    r"""
    Returns the stopping power at energy E [eV] for an openmc material wrt.
    an incident alpha-particle. This is based on the Alsmiller-Estabrook
    correlation [Alsmiller and Estabrook, ORNL-3016, 1960], whence the
    empirical constant 1.866e13 comes.

    Parameters
    ----------
    material: the material through which the alpha-particle is travelling
    E:        the Energy at which the stopping power should be evaluated
    """
    if material is None or E == 0.0:
        return 0
    sumnz = 0
    AE_emp_const = 1.866e13
    for nuc, atoms_per_bcm in material.get_nuclides():
        sumnz += atoms_per_bcm * sqrt(ATOMIC_NUMBER.get(nuc))
    if sumnz != 0:
        return sqrt(E) / (AE_emp_const * sumnz)
    else:
        raise RuntimeError(f"Stopping power of material {material.name} is 0.")
        return 0


def alpha_neutron_yield(material = None, E: float = 0.0) -> float:
    """
    Return the neutron yield per unit volume for an incident alpha particle
    with energy E.

    Parameters
    ----------
    material: the material which the alpha particle is impinging upon
    E:        the energy of the alpha particle
    """
    pass

def alpha_neutron_distribution(material = None, Ea: float = 0.0, dE_E: float = 0.1) -> Tabular:
    """
    Returns a tabular distribution approximating the energy resolved neutron yield from a
    single alpha particle travelling in a material

    Parameters
    ----------
    material: the material which the alpha particle is impinging upon
    Ea:        the energy of the alpha particle
    dE_E:     the relative energy resolution of the resulting distribution
    """
    if material is None:
        return None
    # foreach nuclide in material
    ## get the minimum energy E0, this forms the integration bound together with the emission energy Ea
    ## generate a set of energy points (logscale, with res. dE_E) in [E0,Ea]
    ## foreach energy, E
    ### X-weighted sum over all a,Xn cross sections at E * stopping power(E)
    ## resulting distribution is xs-dist*A.E. xpproximation
    E0=material
    energy=np.logspace(E0,Ea)
