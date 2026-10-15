# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the Brill initial data solver."""

import math

import numpy as np
import pytest
from scipy.special import erf

import brill_initial_data as bid


# --- seed function -------------------------------------------------------


def test_source_is_the_2d_laplacian_of_q():
    params = bid.BrillParams(amplitude=4.7, sigma=1.3)
    rng = np.random.default_rng(1)
    rho, z = rng.uniform(0.1, 3.0, 20), rng.uniform(-3.0, 3.0, 20)
    h = 1e-3
    q = lambda r_, z_: bid.holz_q(params, r_, z_)
    numeric = (
        (q(rho + h, z) - 2 * q(rho, z) + q(rho - h, z))
        + (q(rho, z + h) - 2 * q(rho, z) + q(rho, z - h))
    ) / h ** 2
    assert np.allclose(bid.holz_source(params, rho, z), numeric, rtol=1e-5, atol=1e-7)


def test_source_integrates_to_zero():
    # q vanishes on the axis, so int (q_,rhorho + q_,zz) d^3x = 0.
    grid = bid.Grid(256, 32, 2.0)
    rho, z = grid.mesh()
    source = bid.holz_source(bid.BrillParams(5.0), rho, z)
    w = bid.fejer_half_weights(grid.n_theta)
    integral = grid.h_s * np.sum(grid.dr_ds * grid.r ** 2 * (source @ w))
    assert abs(integral) < 1e-12


# --- quadrature ----------------------------------------------------------


@pytest.mark.parametrize("power", [0, 2, 4, 10])
def test_fejer_weights_integrate_even_polynomials(power):
    n_theta = 16
    x = np.cos(bid.Grid(8, n_theta, 1.0).theta)
    assert np.dot(bid.fejer_half_weights(n_theta), x ** power) == pytest.approx(1.0 / (power + 1), abs=1e-14)


def test_midpoint_even_start():
    h = 1.0 / 64
    s = (np.arange(64) + 0.5) * h
    values = np.cosh(s) + s ** 2
    exact = math.sinh(1.0) + 1.0 / 3.0
    assert bid.midpoint_even_start(values, h) == pytest.approx(exact, abs=1e-11)


# --- operator: manufactured solution ------------------------------------
#
# u = erf(r)/r + r^2 e^{-r^2} P2(cos theta) has the parities of the data
# (even at the origin, monopole 1/r and exponentially small rest at
# infinity) and Laplacian(u) = -(4/sqrt(pi)) e^{-r^2}
#                             + (4 r^4 - 14 r^2) e^{-r^2} P2(cos theta).


def manufactured(grid):
    r = grid.r[:, None]
    p2 = 0.5 * (3.0 * np.cos(grid.theta)[None, :] ** 2 - 1.0)
    gauss = np.exp(-r ** 2)
    u = erf(r) / r + r ** 2 * gauss * p2
    laplacian = -4.0 / math.sqrt(math.pi) * gauss + (4.0 * r ** 4 - 14.0 * r ** 2) * gauss * p2
    return u, laplacian


def manufactured_error(n_s, n_theta):
    grid = bid.Grid(n_s, n_theta, 2.0)
    exact, rhs = manufactured(grid)
    u = bid.solve_linear(grid, np.zeros_like(rhs), rhs)
    return np.max(np.abs(u - exact))


def test_laplacian_converges_at_sixth_order():
    coarse, fine = manufactured_error(64, 16), manufactured_error(128, 32)
    assert fine < 1e-8
    assert math.log2(coarse / fine) > 5.5


def test_flux_mass_of_manufactured_monopole():
    # u -> 1/r means psi = 1 + M/(2r) with M = 2.
    grid = bid.Grid(256, 32, 2.0)
    exact, rhs = manufactured(grid)
    solution = bid.Solution(
        bid.BrillParams(0.0), grid, bid.solve_linear(grid, np.zeros_like(rhs), rhs), np.zeros_like(rhs)
    )
    assert bid.adm_mass(solution).flux == pytest.approx(2.0, abs=1e-8)


# --- Brill data ----------------------------------------------------------


def test_zero_amplitude_is_flat():
    solution = bid.solve_brill(bid.BrillParams(0.0), bid.Grid(64, 16, 2.0))
    assert np.max(np.abs(solution.u)) == 0.0


def test_exponent_factor_one_equals_half_amplitude():
    grid = bid.Grid(64, 16, 2.0)
    u_k1 = bid.solve_brill(bid.BrillParams(5.0, exponent_factor=1), grid).u
    u_k2 = bid.solve_brill(bid.BrillParams(2.5, exponent_factor=2), grid).u
    assert np.allclose(u_k1, u_k2, rtol=0.0, atol=1e-14)


def test_mass_estimators_agree():
    masses = bid.adm_mass(bid.solve_brill(bid.BrillParams(5.0), bid.Grid(256, 32, 2.0)))
    assert masses.flux_spread < 1e-10
    assert masses.source == pytest.approx(masses.flux, abs=1e-8)
    assert masses.brill == pytest.approx(masses.flux, abs=1e-8)


@pytest.mark.parametrize(
    "amplitude, reference",
    [
        (4.6, 0.597),      # Hilditch, Weyhausen, Bruegmann 2017, Table 1
        (4.69, 0.619),
        (4.7, 0.622),
        (4.697, 0.621),
        (6.073, 1.02),     # ibid., Fig. 3
        (5.0, 0.70),       # Khirnov & Ledvinka 2019
        (5.5, 0.84),
    ],
)
def test_adm_mass_matches_published_values(amplitude, reference):
    mass = bid.adm_mass(bid.solve_brill(bid.BrillParams(amplitude), bid.Grid(256, 32, 2.0))).flux
    digits = len(repr(reference).split(".")[1])
    assert round(mass, digits) == pytest.approx(reference, abs=1e-12)


def test_small_amplitude_mass_is_quadratic():
    grid = bid.Grid(128, 16, 2.0)
    ratio = [
        bid.adm_mass(bid.solve_brill(bid.BrillParams(a), grid)).flux / a ** 2 for a in (0.02, 0.01)
    ]
    assert ratio[0] == pytest.approx(ratio[1], rel=2e-3)


# --- output --------------------------------------------------------------


def test_write_table_round_trip(tmp_path):
    h5py = pytest.importorskip("h5py")
    solution = bid.solve_brill(bid.BrillParams(5.0), bid.Grid(64, 16, 2.0))
    masses = bid.adm_mass(solution)
    path = tmp_path / "psi.h5"
    bid.write_table(path, solution, masses)
    with h5py.File(str(path), "r") as handle:
        assert np.array_equal(handle["u"][...], solution.u)
        assert handle.attrs["format"] == "brill-psi-table"
        assert handle.attrs["exponent_factor"] == 2
        assert handle.attrs["n_s"] == 64 and handle.attrs["n_theta"] == 16
        assert handle.attrs["scale_L"] == 2.0
        assert handle.attrs["adm_mass_flux"] == pytest.approx(masses.flux)
