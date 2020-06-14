"""
Ideal gas

:param U: a list of conservative variables [rho, rho u, rho v, rho w, e]
"""

from math import sqrt
from .constant import G, kB_amu


def primitive_to_conservative(rho, v, p):
    """
    Converting primitive variables to a list of conservative variables

    :return U: check the header of this file
    """
    return [
        rho,
        rho * v[0], rho * v[1], rho * v[2],
        .5 * rho * sum([u**2 for u in v]) + p / (G - 1)
    ]


def specific_kinetic_energy(U):
    """
    :return e_kin in (Mpc / Myr)^2
    """
    return .5 * sum([v**2 for v in U[1:4]]) / U[0]**2


def specific_internal_energy(U):
    """
    Internal energy divided by density

    :return e_int in (Mpc / Myr)^2
    """
    return U[4] / U[0] - specific_kinetic_energy(U)


def temperature_over_mu(U):
    """
    :return T / mu in K
    """
    return specific_internal_energy(U) * (G - 1) / kB_amu


def pressure(U):
    """
    :return p in m_p / cm^3 * (Mpc / Myr)^2
    """
    return U[0] * kB_amu * temperature_over_mu(U)


def sound_speed(U):
    """
    :return cs in Mpc / Myr
    """
    return sqrt(G * kB_amu * temperature_over_mu(U))
