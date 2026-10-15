#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Time-symmetric Brill wave initial data: solve the Hamiltonian constraint.

The 3-metric is

    dl^2 = psi^4 [ e^{k q} (drho^2 + dz^2) + rho^2 dphi^2 ],   K_ij = 0,

with k = 2 (the convention of Brill, Holz et al., Hilditch et al.,
Khirnov & Ledvinka and the IDBrillData thorn) or k = 1. The Hamiltonian
constraint reduces to the linear flat-space problem

    Laplacian(psi) + (k/8) (q_,rhorho + q_,zz) psi = 0,   psi -> 1,

where Laplacian is the 3D flat Laplacian in axisymmetry. For the Holz
seed q = A (rho/sigma)^2 exp(-r^2/sigma^2), k = 1 with amplitude A is the
same data as k = 2 with amplitude A/2.

Discretisation
    Polar coordinates (r, theta), compactified as r = L tan(pi s / 2) with
    s in [0, 1). Cell-centred grids in s and in theta in [0, pi/2] keep
    every unknown off the origin, the axis and the equator. The solution
    u = psi - 1 is even in s at s = 0 and even in theta at both ends
    (equatorial and axial symmetry). Outside the Gaussian source it is
    sum_l a_l P_l(cos theta) / r^(l+1) with even l, and 1/r is an odd
    function of 1 - s, so u is odd about s = 1 up to exponentially small
    terms. These parities fill the ghost points of 6th-order centred
    finite differences, and the equation multiplied by r^2,

        a(s) u_ss + b(s) u_s + u_thth + cot(theta) u_th + r^2 V u = r^2 F,

    has coefficients that vanish smoothly at both ends of s. The linear
    system is solved directly with SuperLU.

ADM mass, three estimators
    flux    2 r <u> at large r, where <u> is the monopole (l = 0) part;
    source  (1/2pi) int V psi d^3x, which equals the flux integral when
            the discrete equation holds;
    brill   (1/2pi) int |grad ln psi|^2 d^3x, Brill's positive form,
            which follows from Laplacian(ln psi) = -V - |grad ln psi|^2
            and int q_,rhorho + q_,zz d^3x = 0 (q vanishes on the axis).
    Angular integrals use Fejer's first rule, whose nodes are exactly the
    cell-centred theta grid.

Output: ``write_table`` stores u on the cell-centred (s, theta) grid in an
HDF5 file together with the parameters and the masses, for the
initial-data thorn to interpolate onto the Cactus grid.

Runs with the container's Python 3.8 (numpy, scipy, h5py).
"""

from __future__ import annotations

import argparse
import datetime
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SOLVER_VERSION = "1"
TABLE_FORMAT = "brill-psi-table"
TABLE_FORMAT_VERSION = 1

# 6th-order centred finite-difference weights for offsets -3..3.
D1_WEIGHTS = np.array([-1.0 / 60, 3.0 / 20, -3.0 / 4, 0.0, 3.0 / 4, -3.0 / 20, 1.0 / 60])
D2_WEIGHTS = np.array([1.0 / 90, -3.0 / 20, 3.0 / 2, -49.0 / 18, 3.0 / 2, -3.0 / 20, 1.0 / 90])
OFFSETS = np.arange(-3, 4)


# ---------------------------------------------------------------------
# Seed function
# ---------------------------------------------------------------------


@dataclass
class BrillParams:
    """Parameters of the Holz-type Brill wave.

    Attributes:
        amplitude: A in q = A (rho/sigma)^2 exp(-r^2/sigma^2).
        sigma: Width sigma.
        exponent_factor: k in e^{k q}; 2 is the standard convention.
    """

    amplitude: float
    sigma: float = 1.0
    exponent_factor: int = 2

    def __post_init__(self) -> None:
        if self.exponent_factor not in (1, 2):
            raise ValueError("exponent_factor must be 1 or 2")
        if self.sigma <= 0.0:
            raise ValueError("sigma must be positive")


def holz_q(params: BrillParams, rho: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Holz seed q = A (rho/sigma)^2 exp(-(rho^2 + z^2)/sigma^2)."""
    s2 = params.sigma ** 2
    return params.amplitude * rho ** 2 / s2 * np.exp(-(rho ** 2 + z ** 2) / s2)


def holz_source(params: BrillParams, rho: np.ndarray, z: np.ndarray) -> np.ndarray:
    """q_,rhorho + q_,zz for the Holz seed.

    Equals (2A/sigma^2) e^{-r^2/sigma^2} [1 - 6 rho^2/sigma^2 + 2 rho^2 r^2/sigma^4].
    """
    s2 = params.sigma ** 2
    r2 = rho ** 2 + z ** 2
    return (
        2.0 * params.amplitude / s2 * np.exp(-r2 / s2)
        * (1.0 - 6.0 * rho ** 2 / s2 + 2.0 * rho ** 2 * r2 / s2 ** 2)
    )


# ---------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------


@dataclass
class Grid:
    """Cell-centred compactified polar grid.

    Attributes:
        n_s: Cells in s over [0, 1).
        n_theta: Cells in theta over [0, pi/2].
        scale: L in r = L tan(pi s / 2).
    """

    n_s: int
    n_theta: int
    scale: float
    s: np.ndarray = field(init=False, repr=False)
    theta: np.ndarray = field(init=False, repr=False)
    h_s: float = field(init=False)
    h_theta: float = field(init=False)

    def __post_init__(self) -> None:
        if self.n_s < 8 or self.n_theta < 4:
            raise ValueError("grid too small for 6th-order stencils")
        self.h_s = 1.0 / self.n_s
        self.h_theta = 0.5 * math.pi / self.n_theta
        self.s = (np.arange(self.n_s) + 0.5) * self.h_s
        self.theta = (np.arange(self.n_theta) + 0.5) * self.h_theta

    @property
    def r(self) -> np.ndarray:
        """Radius at the s cell centres."""
        return radius_of_s(self.s, self.scale)

    @property
    def dr_ds(self) -> np.ndarray:
        """dr/ds at the s cell centres."""
        phi = 0.5 * math.pi * self.s
        return 0.5 * math.pi * self.scale / np.cos(phi) ** 2

    def mesh(self) -> Tuple[np.ndarray, np.ndarray]:
        """(rho, z) on the (n_s, n_theta) grid."""
        r = self.r[:, None]
        return r * np.sin(self.theta)[None, :], r * np.cos(self.theta)[None, :]


def radius_of_s(s: np.ndarray, scale: float) -> np.ndarray:
    """r = L tan(pi s / 2)."""
    return scale * np.tan(0.5 * math.pi * np.asarray(s))


def s_of_radius(r: np.ndarray, scale: float) -> np.ndarray:
    """Inverse of ``radius_of_s``: s = (2/pi) arctan(r / L)."""
    return 2.0 / math.pi * np.arctan(np.asarray(r) / scale)


def radial_coefficients(grid: Grid) -> Tuple[np.ndarray, np.ndarray]:
    """Coefficients of r^2 (u_rr + 2 u_r / r) = a u_ss + b u_s.

    For r = L tan(phi), phi = pi s / 2:
    a = (4/pi^2) sin^2(phi) cos^2(phi), b = (4/pi) sin(phi) cos^3(phi).
    """
    phi = 0.5 * math.pi * grid.s
    sin, cos = np.sin(phi), np.cos(phi)
    return 4.0 / math.pi ** 2 * sin ** 2 * cos ** 2, 4.0 / math.pi * sin * cos ** 3


def fejer_half_weights(n_theta: int) -> np.ndarray:
    """Weights w_j with int_0^1 g(x) dx ~ sum_j w_j g(cos theta_j).

    Fejer's first rule on [-1, 1] with 2 n_theta Chebyshev nodes
    theta_j = (j + 1/2) pi / (2 n_theta); the nodes in [0, pi/2] are the
    cell-centred theta grid. Valid for g symmetric about x = 0.
    """
    n = 2 * n_theta
    theta = (np.arange(n_theta) + 0.5) * math.pi / n
    k = np.arange(1, n // 2 + 1)
    series = np.cos(2.0 * np.outer(theta, k)) / (4.0 * k ** 2 - 1.0)
    return 2.0 / n * (1.0 - 2.0 * series.sum(axis=1))


# ---------------------------------------------------------------------
# Operator and solver
# ---------------------------------------------------------------------


def _mirror(index: np.ndarray, n: int, upper_parity: float) -> Tuple[np.ndarray, np.ndarray]:
    """Map stencil indices outside [0, n) onto the grid with their parities.

    Cell-centred mirror at the lower face (even) and at the upper face
    (``upper_parity``).
    """
    mapped = index.copy()
    parity = np.ones_like(index, dtype=float)
    low = index < 0
    mapped[low] = -index[low] - 1
    high = index >= n
    mapped[high] = 2 * n - 1 - index[high]
    parity[high] = upper_parity
    return mapped, parity


def build_operator(grid: Grid, potential: np.ndarray) -> sp.csr_matrix:
    """Sparse matrix of r^2 (Laplacian + V) acting on u.

    Args:
        grid: The grid.
        potential: V on the (n_s, n_theta) grid.

    Returns:
        The (n_s n_theta)^2 operator, rows and columns in C order of (i, j).
    """
    n_s, n_t = grid.n_s, grid.n_theta
    a, b = radial_coefficients(grid)
    cot = 1.0 / np.tan(grid.theta)
    r2 = grid.r ** 2
    rows: List[np.ndarray] = []
    cols: List[np.ndarray] = []
    vals: List[np.ndarray] = []
    ii, jj = np.meshgrid(np.arange(n_s), np.arange(n_t), indexing="ij")
    row = (ii * n_t + jj).ravel()

    for k, w1, w2 in zip(OFFSETS, D1_WEIGHTS, D2_WEIGHTS):
        # s direction: odd parity about s = 1.
        target, parity = _mirror(ii + k, n_s, -1.0)
        coeff = (a[:, None] * w2 / grid.h_s ** 2 + b[:, None] * w1 / grid.h_s) * np.ones((1, n_t))
        rows.append(row)
        cols.append((target * n_t + jj).ravel())
        vals.append((coeff * parity).ravel())
        # theta direction: even parity at both ends.
        target, parity = _mirror(jj + k, n_t, 1.0)
        coeff = (w2 / grid.h_theta ** 2 + cot[None, :] * w1 / grid.h_theta) * np.ones((n_s, 1))
        rows.append(row)
        cols.append((ii * n_t + target).ravel())
        vals.append((coeff * parity).ravel())

    rows.append(row)
    cols.append(row)
    vals.append((r2[:, None] * potential).ravel())
    size = n_s * n_t
    return sp.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(size, size)
    )


def solve_linear(grid: Grid, potential: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve (Laplacian + V) u = F with u -> 0 at infinity.

    Args:
        grid: The grid.
        potential: V on the grid.
        rhs: F on the grid.

    Returns:
        u on the (n_s, n_theta) grid.
    """
    operator = build_operator(grid, potential).tocsc()
    rhs_scaled = (grid.r[:, None] ** 2 * rhs).ravel()
    return spla.spsolve(operator, rhs_scaled).reshape(grid.n_s, grid.n_theta)


@dataclass
class Solution:
    """Solved Brill initial data.

    Attributes:
        params: Wave parameters.
        grid: Grid.
        u: psi - 1 on the cell-centred grid, shape (n_s, n_theta).
        potential: V = (k/8)(q_,rhorho + q_,zz) on the grid.
    """

    params: BrillParams
    grid: Grid
    u: np.ndarray
    potential: np.ndarray

    @property
    def psi(self) -> np.ndarray:
        return 1.0 + self.u


def solve_brill(params: BrillParams, grid: Grid) -> Solution:
    """Solve the Hamiltonian constraint for one Brill wave.

    Args:
        params: Wave parameters.
        grid: Grid.

    Returns:
        The solution.
    """
    rho, z = grid.mesh()
    potential = params.exponent_factor / 8.0 * holz_source(params, rho, z)
    u = solve_linear(grid, potential, -potential)
    return Solution(params, grid, u, potential)


# ---------------------------------------------------------------------
# Derivatives, quadrature and the ADM mass
# ---------------------------------------------------------------------


def derivative(field_: np.ndarray, axis: int, h: float, upper_parity: float, order: int) -> np.ndarray:
    """6th-order centred derivative of a grid function, ghosts by parity.

    Args:
        field_: Values on the (n_s, n_theta) grid.
        axis: 0 for s, 1 for theta.
        h: Grid spacing along ``axis``.
        upper_parity: Parity at the upper face of ``axis``.
        order: 1 or 2.

    Returns:
        The derivative on the same grid.
    """
    weights = D1_WEIGHTS if order == 1 else D2_WEIGHTS
    n = field_.shape[axis]
    result = np.zeros_like(field_)
    base = np.arange(n)
    for k, w in zip(OFFSETS, weights):
        if w == 0.0:
            continue
        target, parity = _mirror(base + k, n, upper_parity)
        shifted = np.take(field_, target, axis=axis)
        shape = [1, 1]
        shape[axis] = n
        result += w * parity.reshape(shape) * shifted
    return result / h ** order


def midpoint_even_start(values: np.ndarray, h: float, fit_points: int = 10, degree: int = 7) -> float:
    """int_0^1 F(s) ds from cell-centred samples, for F even about s = 0.

    The midpoint rule is spectrally accurate at an even endpoint; at
    s = 1 the Euler-Maclaurin terms are added with F'(1) and F'''(1) from
    a polynomial fit to the last cells.
    """
    n = values.size
    s = (np.arange(n) + 0.5) * h
    total = h * values.sum()
    coeffs = np.polyfit(s[-fit_points:] - 1.0, values[-fit_points:], degree)
    poly = np.poly1d(coeffs)
    d1 = poly.deriv(1)(0.0)
    d3 = poly.deriv(3)(0.0)
    return total + h ** 2 / 24.0 * d1 - 7.0 * h ** 4 / 5760.0 * d3


@dataclass
class MassEstimates:
    """ADM mass from three independent routes.

    Attributes:
        flux: 2 r <u> averaged over the far zone.
        flux_spread: Max minus min of 2 r <u> over the far zone.
        source: (1/2pi) int V psi d^3x.
        brill: (1/2pi) int |grad ln psi|^2 d^3x.
    """

    flux: float
    flux_spread: float
    source: float
    brill: float


def adm_mass(solution: Solution, far_zone: Tuple[float, float] = (20.0, 200.0)) -> MassEstimates:
    """Compute the three ADM mass estimates.

    Args:
        solution: Solved data.
        far_zone: Radius range for the flux estimate (outside the source).

    Returns:
        The estimates.
    """
    grid = solution.grid
    w = fejer_half_weights(grid.n_theta)
    r, dr_ds = grid.r, grid.dr_ds
    monopole = solution.u @ w  # int_0^1 u dx: the l = 0 coefficient
    zone = (r >= far_zone[0]) & (r <= far_zone[1])
    if not zone.any():
        raise ValueError("no grid points in the far zone; refine the grid or widen the zone")
    flux_values = 2.0 * r[zone] * monopole[zone]

    source = 2.0 * grid.h_s * float(np.sum(dr_ds * r ** 2 * ((solution.potential * solution.psi) @ w)))

    psi = solution.psi
    psi_s = derivative(solution.u, 0, grid.h_s, -1.0, 1)
    psi_t = derivative(solution.u, 1, grid.h_theta, 1.0, 1)
    integrand = (psi_s ** 2 * (r ** 2 / dr_ds)[:, None] + psi_t ** 2 * dr_ds[:, None]) / psi ** 2
    brill = 2.0 * midpoint_even_start(integrand @ w, grid.h_s)
    return MassEstimates(
        flux=float(np.mean(flux_values)),
        flux_spread=float(np.ptp(flux_values)),
        source=source,
        brill=brill,
    )


# ---------------------------------------------------------------------
# Table interpolation (mirrors thorns/BrillCollapse/BrillPsiTable)
# ---------------------------------------------------------------------


def _lagrange_weights(x: np.ndarray, first: np.ndarray, points: int) -> np.ndarray:
    """Lagrange weights at x for the integer nodes first..first+points-1."""
    nodes = first[:, None] + np.arange(points)[None, :]
    weights = np.ones((x.size, points))
    for k in range(points):
        for m in range(points):
            if m != k:
                weights[:, k] *= (x - nodes[:, m]) / float(k - m)
    return weights


def evaluate_table(
    u: np.ndarray,
    scale: float,
    r: np.ndarray,
    theta: np.ndarray,
    points: int = 8,
    parity_s_upper: float = -1.0,
) -> np.ndarray:
    """Interpolate a u table at (r, theta) exactly as the Cactus thorn does.

    Tensor-product Lagrange interpolation in (s, theta) with ``points``
    nodes per direction centred on the target, table ghosts filled by
    parity (even at the lower faces and at theta = pi/2, ``parity_s_upper``
    at s = 1).

    Args:
        u: Table on the cell-centred grid, shape (n_s, n_theta).
        scale: L in r = L tan(pi s / 2).
        r: Radii.
        theta: Polar angles in [0, pi/2].
        points: Stencil width per direction (even).
        parity_s_upper: Parity of u about s = 1.

    Returns:
        u at the points, shaped like ``r``.
    """
    n_s, n_t = u.shape
    r_flat = np.asarray(r, dtype=float).ravel()
    t_flat = np.asarray(theta, dtype=float).ravel()
    xs = s_of_radius(r_flat, scale) * n_s - 0.5
    xt = t_flat / (0.5 * math.pi / n_t) - 0.5
    half = points // 2
    fs = np.floor(xs).astype(int) - half + 1
    ft = np.floor(xt).astype(int) - half + 1
    ws = _lagrange_weights(xs, fs, points)
    wt = _lagrange_weights(xt, ft, points)
    offsets = np.arange(points)
    i, ps = _mirror(fs[:, None] + offsets[None, :], n_s, parity_s_upper)
    j, pt = _mirror(ft[:, None] + offsets[None, :], n_t, 1.0)
    values = u[i[:, :, None], j[:, None, :]] * ps[:, :, None] * pt[:, None, :]
    return np.einsum("pa,pb,pab->p", ws, wt, values).reshape(np.shape(r))


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------


def write_table(path: Path, solution: Solution, masses: MassEstimates) -> None:
    """Write u on the (s, theta) grid with its metadata to HDF5.

    Layout: dataset ``u`` of shape (n_s, n_theta) at s_i = (i + 1/2)/n_s,
    theta_j = (j + 1/2) (pi/2)/n_theta, with r = L tan(pi s/2) and theta
    measured from the +z axis. Ghost parities: even at s = 0, odd at
    s = 1, even at theta = 0 and pi/2.

    Args:
        path: Output file.
        solution: Solved data.
        masses: ADM mass estimates.
    """
    import h5py  # deferred: only the writer needs it

    grid, params = solution.grid, solution.params
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(str(path), "w") as handle:
        dataset = handle.create_dataset("u", data=solution.u)
        # Strings are fixed-length ASCII so that the C reader can use them
        # without variable-length string handling.
        dataset.attrs["description"] = np.bytes_("psi - 1 on the cell-centred compactified polar grid")
        attrs = handle.attrs
        attrs["format"] = np.bytes_(TABLE_FORMAT)
        attrs["format_version"] = TABLE_FORMAT_VERSION
        attrs["solver_version"] = np.bytes_(SOLVER_VERSION)
        attrs["created_at"] = np.bytes_(
            datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        )
        attrs["family"] = np.bytes_("holz")
        attrs["amplitude"] = params.amplitude
        attrs["sigma"] = params.sigma
        attrs["exponent_factor"] = params.exponent_factor
        attrs["mapping"] = np.bytes_("r = L tan(pi s / 2)")
        attrs["scale_L"] = grid.scale
        attrs["n_s"] = grid.n_s
        attrs["n_theta"] = grid.n_theta
        attrs["parity_s_lower"] = 1
        attrs["parity_s_upper"] = -1
        attrs["parity_theta_lower"] = 1
        attrs["parity_theta_upper"] = 1
        attrs["adm_mass_flux"] = masses.flux
        attrs["adm_mass_source"] = masses.source
        attrs["adm_mass_brill"] = masses.brill
        attrs["psi_min"] = float(solution.psi.min())


# ---------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------


def report(params: BrillParams, grid: Grid) -> Tuple[Solution, MassEstimates]:
    """Solve and print one line of diagnostics."""
    solution = solve_brill(params, grid)
    masses = adm_mass(solution)
    print(
        "A={:<8g} k={} n_s={:<5d} n_theta={:<4d} L={:<4g} "
        "M_flux={:.12f} (spread {:.1e})  M_source={:.12f}  M_brill={:.12f}  psi_min={:.6f}".format(
            params.amplitude, params.exponent_factor, grid.n_s, grid.n_theta, grid.scale,
            masses.flux, masses.flux_spread, masses.source, masses.brill, float(solution.psi.min()),
        )
    )
    return solution, masses


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--amplitude", "-A", type=float, required=True, help="Brill amplitude A")
    parser.add_argument("--sigma", type=float, default=1.0, help="width sigma (default: 1)")
    parser.add_argument(
        "--exponent-factor", type=int, choices=(1, 2), default=2,
        help="k in e^{k q} (default: 2, the standard convention)",
    )
    parser.add_argument("--n-s", type=int, default=512, help="cells in s (default: 512)")
    parser.add_argument("--n-theta", type=int, default=64, help="cells in theta (default: 64)")
    parser.add_argument("--scale", "-L", type=float, default=2.0, help="L in r = L tan(pi s/2) (default: 2)")
    parser.add_argument(
        "--convergence", action="store_true",
        help="solve at n/4, n/2 and n and report the convergence of the masses",
    )
    parser.add_argument("--out", type=Path, help="write the psi table to this HDF5 file")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    params = BrillParams(args.amplitude, args.sigma, args.exponent_factor)
    if args.convergence:
        results = []
        for factor in (4, 2, 1):
            grid = Grid(args.n_s // factor, args.n_theta // factor, args.scale)
            results.append(report(params, grid)[1])
        for name in ("flux", "source", "brill"):
            m = [getattr(r, name) for r in results]
            d1, d2 = m[1] - m[0], m[2] - m[1]
            order = math.log2(abs(d1 / d2)) if d2 != 0.0 and d1 != 0.0 else float("nan")
            print("  {:<6}: differences {:.2e}, {:.2e}  observed order {:.2f}".format(name, d1, d2, order))
        return 0
    solution, masses = report(params, Grid(args.n_s, args.n_theta, args.scale))
    if args.out is not None:
        write_table(args.out, solution, masses)
        print("wrote {}".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
