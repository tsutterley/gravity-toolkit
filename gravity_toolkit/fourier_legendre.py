#!/usr/bin/env python
"""
fourier_legendre.py
Original IDL code gen_plms.pro written by Sean Swenson
Adapted by Tyler Sutterley (07/2026)

Computes Fourier coefficients of the associated Legendre functions

CALLING SEQUENCE:
    Almk = fourier_legendre(lmax,mmax)

INPUTS:
    lmax: maximum spherical harmonic degree
    mmax: maximum spherical harmonic order

OUTPUTS:
    Almk: Fourier coefficients

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 07/2027: add citations to docstrings
    Updated 03/2023: improve typing for variables in docstrings
    Updated 10/2022: add polynomial function for calculating gradients
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: cleaned up program for public release
    Updated 07/2020: added function docstrings
    Updated 06/2019: using Python3 compatible division
    Written 04/2013
"""

from __future__ import division
import numpy as np


def fourier_legendre(lmax, mmax):
    """
    Computes Fourier coefficients of the associated Legendre functions
    :cite:p:`Hofsommer:1960wg,Gruber:2016hn`

    Parameters
    ----------
    lmax: int
        Upper bound of Spherical Harmonic Degrees
    mmax: int
        Upper bound of Spherical Harmonic Orders

    Returns
    -------
    Almk: np.ndarray
        Fourier coefficients
    """

    # allocate for output fourier coefficients
    Almk = np.zeros((lmax + 1, lmax + 1, lmax + 1))
    l_even = np.arange(0, lmax + 1, 2)
    l_odd = np.arange(1, lmax, 2)
    m_even = np.arange(0, mmax + 1, 2)
    m_odd = np.arange(1, mmax, 2)

    # First compute m=0, m=1 terms
    # Compute m = 0, l = even terms
    Almk[l_even, 0, 0] = 1.0
    a1 = (l_even * (l_even + 1.0)) * Almk[l_even, 0, 0]
    Almk[l_even, 0, 2] = a1 / (l_even * (l_even + 1.0) - 2.0)
    for j in range(2, lmax, 2):  # equivalent to 2:lmax-2
        a1 = 2.0 * (l_even * (l_even + 1.0) - j**2.0) * Almk[l_even, 0, j]
        a2 = ((j - 2.0) * (j - 1.0) - l_even * (l_even + 1.0)) * Almk[
            l_even, 0, j - 2
        ]
        dfactor = l_even * (l_even + 1.0) - (j + 2.0) * (j + 1.0)
        Almk[l_even, 0, j + 2] = (a1 + a2) / dfactor

    # Special case for j = 0 fourier coefficient
    Almk[l_even, 0, 0] = Almk[l_even, 0, 0] / 2.0
    # Normalize overall sum to 2 for m == 0
    norm = np.zeros((len(l_even)))
    for j in range(0, lmax + 2, 2):  # equivalent to 0:lmax
        ptemp = np.squeeze(Almk[l_even[:, np.newaxis], 0, m_even])
        dtemp = (
            1.0 / (1.0 - j - m_even)
            + 1.0 / (1.0 + j - m_even)
            + 1.0 / (1.0 - j + m_even)
            + 1.0 / (1.0 + j + m_even)
        )
        norm[l_even // 2] = (
            norm[l_even // 2] + Almk[l_even, 0, j] * np.dot(ptemp, dtemp) / 2.0
        )
    # normalize Almks
    norm = np.sqrt(norm / 2.0)
    for l in range(0, lmax + 2, 2):  # equivalent to 0:lmax
        Almk[l, 0, :] = Almk[l, 0, :] / norm[l // 2]

    # Compute m = 0, l = odd terms
    Almk[l_odd, 0, 1] = 1.0
    a1 = (2.0 - l_odd * (l_odd + 1.0)) * Almk[l_odd, 0, 1]
    Almk[l_odd, 0, 3] = a1 / (6.0 - l_odd * (l_odd + 1.0))
    for j in range(3, lmax - 1, 2):  # equivalent to 3:lmax-3
        a1 = 2.0 * (l_odd * (l_odd + 1.0) - j**2.0) * Almk[l_odd, 0, j]
        a2 = ((j - 2.0) * (j - 1.0) - l_odd * (l_odd + 1.0)) * Almk[
            l_odd, 0, j - 2
        ]
        dfactor = l_odd * (l_odd + 1.0) - (j + 2.0) * (j + 1.0)
        Almk[l_odd, 0, j + 2] = (a1 + a2) / dfactor

    # Normalize overall sum to 2 for m == 0
    norm = np.zeros((len(l_odd)))
    for j in range(1, lmax + 1, 2):  # equivalent to 1:lmax-1
        ptemp = np.squeeze(Almk[l_odd[:, np.newaxis], 0, m_odd])
        dtemp = (
            1.0 / (1.0 - j - m_odd)
            + 1.0 / (1.0 + j - m_odd)
            + 1.0 / (1.0 - j + m_odd)
            + 1.0 / (1.0 + j + m_odd)
        )
        norm[(l_odd - 1) // 2] = (
            norm[(l_odd - 1) // 2]
            + Almk[l_odd, 0, j] * np.dot(ptemp, dtemp) / 2.0
        )
    # normalize Almks
    norm = np.sqrt(norm / 2.0)
    for l in range(1, lmax + 1, 2):  # equivalent to 1:lmax-1
        Almk[l, 0, :] = Almk[l, 0, :] / norm[(l - 1) // 2]

    # Compute m = 1, l = even terms
    Almk[l_even, 1, 0] = 0.0
    Almk[l_even, 1, 2] = 1.0
    for j in range(2, lmax, 2):  # equivalent to 2:lmax-2
        a1 = 2.0 * (l_even * (l_even + 1) - j**2.0 - 2.0) * Almk[l_even, 1, j]
        a2 = ((j - 2.0) * (j - 1.0) - l_even * (l_even + 1)) * Almk[
            l_even, 1, j - 2
        ]
        dfactor = l_even * (l_even + 1.0) - (j + 2.0) * (j + 1.0)
        Almk[l_even, 1, j + 2] = (a1 + a2) / dfactor

    # Normalize overall sum to 4 for m == 1
    # different norm than that of the cosine series
    norm = np.zeros((len(l_even)))
    for j in range(0, lmax + 2, 2):  # equivalent to 0:lmax
        ptemp = np.squeeze(Almk[l_even[:, np.newaxis], 1, m_even])
        dtemp = (
            -1.0 / (1.0 - j - m_even)
            + 1.0 / (1 + j - m_even)
            + 1.0 / (1.0 - j + m_even)
            - 1.0 / (1 + j + m_even)
        )
        norm[l_even // 2] = (
            norm[l_even // 2] + Almk[l_even, 1, j] * np.dot(ptemp, dtemp) / 2.0
        )
    # normalize Almks
    norm = np.sqrt(norm / 4.0)
    for l in range(0, lmax + 2, 2):  # equivalent to 0:lmax
        Almk[l, 1, :] = Almk[l, 1, :] / norm[l // 2]

    # Compute m = 1, l = odd terms
    Almk[l_odd, 1, 1] = 1.0
    Almk[l_odd, 1, 3] = (
        3.0
        * (l_odd * (l_odd + 1) - 2)
        * Almk[l_odd, 1, 1]
        / (l_odd * (l_odd + 1) - 6)
    )
    for j in range(3, lmax - 1, 2):  # equivalent to 3:lmax-3
        a1 = 2.0 * (l_odd * (l_odd + 1.0) - j**2.0 - 2.0) * Almk[l_odd, 1, j]
        a2 = ((j - 2.0) * (j - 1.0) - l_odd * (l_odd + 1.0)) * Almk[
            l_odd, 1, j - 2
        ]
        dfactor = l_odd * (l_odd + 1.0) - (j + 2.0) * (j + 1.0)
        Almk[l_odd, 1, j + 2] = (a1 + a2) / dfactor

    # Normalize overall sum to 4 for m == 1
    norm = np.zeros((len(l_odd)))
    for j in range(1, lmax + 1, 2):  # equivalent to 1:lmax-1
        ptemp = np.squeeze(Almk[l_odd[:, np.newaxis], 1, m_odd])
        dtemp = (
            -1.0 / (1.0 - j - m_odd)
            + 1.0 / (1.0 + j - m_odd)
            + 1.0 / (1.0 - j + m_odd)
            - 1.0 / (1.0 + j + m_odd)
        )
        norm[(l_odd - 1) // 2] = (
            norm[(l_odd - 1) // 2]
            + Almk[l_odd, 1, j] * np.dot(ptemp, dtemp) / 2.0
        )
    # normalize Almks
    norm = np.sqrt(norm / 4.0)
    for l in range(1, lmax + 1, 2):  # equivalent to 1:lmax-1
        Almk[l, 1, :] = Almk[l, 1, :] / norm[(l - 1) // 2]

    # Compute coefficients for m > 0
    # m = 0 terms on rhs have different normalization
    m = 0
    # m = 0, l = even terms
    for l in range(m, lmax - 1):  # equivalent to m:lmax-2
        a1 = (
            np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
            * Almk[l, m, m_even]
        )
        a2 = (
            np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
            * Almk[l + 2, m, m_even]
        )
        a3 = (
            np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0) / 2.0)
            * Almk[l, m + 2, m_even]
        )
        dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0) / 2.0)
        Almk[l + 2, m + 2, m_even] = (a1 - a2 + a3) / dfactor

    # m = 0, l = odd terms
    for l in range(m + 1, lmax - 1):  # equivalent to m+1:lmax-2
        a1 = (
            np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
            * Almk[l, m, m_odd]
        )
        a2 = (
            np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
            * Almk[l + 2, m, m_odd]
        )
        a3 = (
            np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0) / 2.0)
            * Almk[l, m + 2, m_odd]
        )
        dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0) / 2.0)
        Almk[l + 2, m + 2, m_odd] = (a1 - a2 + a3) / dfactor

    # m = even terms
    for m in range(2, lmax, 2):  # equivalent to 2:lmax-2
        # m = even, > 2, l = even terms
        for l in range(m, lmax, 2):  # equivalent to m:lmax-2
            a1 = (
                np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
                * Almk[l, m, m_even]
            )
            a2 = (
                np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
                * Almk[l + 2, m, m_even]
            )
            a3 = (
                np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0))
                * Almk[l, m + 2, m_even]
            )
            dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0))
            Almk[l + 2, m + 2, m_even] = (a1 - a2 + a3) / dfactor

        # m = even, > 2, l = odd terms
        for l in range(m + 1, lmax - 1, 2):
            a1 = (
                np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
                * Almk[l, m, m_odd]
            )
            a2 = (
                np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
                * Almk[l + 2, m, m_odd]
            )
            a3 = (
                np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0))
                * Almk[l, m + 2, m_odd]
            )
            dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0))
            Almk[l + 2, m + 2, m_odd] = (a1 - a2 + a3) / dfactor

    # m = odd terms
    for m in range(1, lmax - 1, 2):  # equivalent to 1:lmax-3
        # m = odd, > 1, l = even terms
        for l in range(m + 1, lmax - 1, 2):  # equivalent to m+1,lmax-2
            a1 = (
                np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
                * Almk[l, m, m_even]
            )
            a2 = (
                np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
                * Almk[l + 2, m, m_even]
            )
            a3 = (
                np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0))
                * Almk[l, m + 2, m_even]
            )
            dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0))
            Almk[l + 2, m + 2, m_even] = (a1 - a2 + a3) / dfactor

        # m = odd, > 1, l = odd terms
        for l in range(m, lmax - 1, 2):  # equivalent to m:lmax-2
            a1 = (
                np.sqrt((l + m + 2.0) * (l + m + 1.0) / (2.0 * l + 1.0))
                * Almk[l, m, m_odd]
            )
            a2 = (
                np.sqrt((l - m + 1.0) * (l - m + 2.0) / (2.0 * l + 5.0))
                * Almk[l + 2, m, m_odd]
            )
            a3 = (
                np.sqrt((l - m) * (l - m - 1.0) / (2.0 * l + 1.0))
                * Almk[l, m + 2, m_odd]
            )
            dfactor = np.sqrt((l + m + 4.0) * (l + m + 3.0) / (2.0 * l + 5.0))
            Almk[l + 2, m + 2, m_odd] = (a1 - a2 + a3) / dfactor

    # return the fourier coefficients
    return Almk


def legendre_gradient(lmax, mmax):
    """
    Calculates functions for evaluating the integral of a
    product of two fourier series

    Parameters
    ----------
    lmax: int
        Upper bound of Spherical Harmonic Degrees
    mmax: int
        Upper bound of Spherical Harmonic Orders

    Returns
    -------
    Vlmk: np.ndarray
        Fourier coefficients for meridional gradients
    Wlmk: np.ndarray
        Fourier coefficients for zonal gradients
    """
    # compute the fourier coefficients of the associated legendre functions
    Almk = fourier_legendre(lmax, mmax)
    # allocate for output fourier coefficients
    Vlmk = np.zeros((lmax + 1, lmax + 1, lmax + 1))
    Wlmk = np.zeros((lmax + 1, lmax + 1, lmax + 1))
    # for each spherical harmonic degree
    for l in range(1, lmax + 1):
        # degree dependent factor
        dfactor = np.sqrt((2.0 * l + 1.0) / (2.0 * l - 1.0))
        # m=0 special case
        Vfact = np.sqrt(l * (l + 1.0) / 2.0)
        Vlmk[l, 0, :] = Vfact * Almk[l, 1, :]
        for m in range(2, l + 1):  # from 2 to l
            Vfact = np.sqrt((l + m) * (l - m + 1.0) / 4.0)
            Wfact = dfactor * np.sqrt((l - m) * (l - m + 1) / 4.0)
            Vlmk[l, m - 1, :] = Vfact * Almk[l, m, :]
            Wlmk[l, m - 1, :] = -Wfact * Almk[l - 1, m, :]
        # m = 1 terms
        Vfact = np.sqrt(l * (l + 1.0) / 2.0)
        Wfact = dfactor * np.sqrt(l * (l + 1.0) / 2.0)
        Vlmk[l, 1, :] -= Vfact * Almk[l, 0, :]
        Wlmk[l, 1, :] += dfactor * Wfact * Almk[l - 1, 0, :]
        for m in range(2, l + 1):  # from 2 to l
            Vfact = np.sqrt((l + m) * (l - m + 1.0) / 4.0)
            Wfact = dfactor * np.sqrt((l + m) * (l + m - 1) / 4.0)
            Vlmk[l, m, :] -= Vfact * Almk[l, m - 1, :]
            Wlmk[l, m, :] += Wfact * Almk[l - 1, m - 1, :]
        # normalizations
        Vlmk[l, :, :] /= np.sqrt(l * (l + 1.0))
        Wlmk[l, :, :] /= np.sqrt(l * (l + 1.0))
    # return the coefficients
    return (Vlmk, Wlmk)
