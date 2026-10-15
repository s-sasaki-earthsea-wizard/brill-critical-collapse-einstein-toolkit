# Brill Wave Vacuum Critical Collapse with the Einstein Toolkit

Numerical study of the threshold of black hole formation from pure
gravitational waves: time-symmetric Brill waves evolved with the
Einstein Toolkit (ET). A bisection in the wave amplitude `A` locates the
critical amplitude `A*`, and the apparent horizon masses `M_AH` of the
supercritical runs give the critical exponent `gamma` through the
scaling law

```text
M_AH ∝ (A - A*)^gamma
```

> **Status: skeleton.** The environment, the parfile template and the
> script interfaces are in place, but the project cannot produce physics
> yet: the initial data solver is an open design decision (see
> [Open issues](#open-issues)), and the bisection, monitoring and
> analysis logic is still to be implemented.

## Physical setup

| Item | Choice |
| --- | --- |
| Initial data | Time-symmetric Brill wave, Holz form `q = A rho^2 exp(-r^2 / sigma^2)`, `sigma = 1` |
| Search range | `A ∈ [4.0, 6.0]` |
| Evolution | McLachlan `ML_BSSN`, RK4, 4th-order finite differences |
| Gauge | 1+log slicing with lapse advection, Gamma-driver shift |
| Symmetry | Octant (`x, y, z >= 0`) via `ReflectionSymmetry` |
| Grid | Carpet, 4 fixed refinement levels centred on the origin, finest `dx = 0.025`, outer boundary at 24 |
| Dissipation | Kreiss-Oliger, 5th order (`Dissipation` thorn) |
| Horizons | `AHFinderDirect` every 16 iterations, initial guess radius 0.2 |
| End time | `t_max = 30` |

## Prerequisites

- Linux host with Docker and Docker Compose v2.
- The GW230529 project's Einstein Toolkit image (`gw230529-et:local`),
  which this project overlays. Build it once in
  [gw230529-einstein-toolkit](https://github.com/s-sasaki-earthsea-wizard/gw230529-einstein-toolkit):

  ```bash
  cd ../gw230529-einstein-toolkit && make docker-build
  ```

  It already contains every thorn used here (ET_2025_05 "Kruskal",
  MPICH only); `make docker-check` verifies that.
- About 12 GB of RAM per run for the default grid (Carpet's own
  estimate is 10.9 GB) and a local SSD for the output directory.

## Quick start

```bash
make docker-setup        # create .env and the output directory
make docker-build        # overlay build on top of gw230529-et:local
make docker-up           # start the container
make docker-check        # thorns, OMP_NUM_THREADS=1, MPICH linkage, MPI communicator
make check               # byte-compile the scripts with the container's Python 3.8
make test                # pytest suite
make bisection-dry-run   # render the first parfile and print the run plan
```

Once the initial data and the runner are implemented:

```bash
make single-run AMP=4.0               # smoke test on the subcritical side
make bisection                        # full campaign (A_MIN, A_MAX, MAX_ITER, CAMPAIGN)
make status                           # current amplitude, ETA, memory
make analyze                          # gamma and reports/critical_scaling.png
```

## Workflow

1. **Bisection** (`scripts/bisection_runner.py`). For each step the
   runner renders the template at the bracket midpoint, runs
   `cactus_sim` under `mpirun.mpich`, and classifies the run:
   - `collapse`: `AHFinderDirect` found a horizon with `M_AH > 0`, so the
     upper bound moves down;
   - `dispersion`: no horizon, and the central lapse recovered, so the
     lower bound moves up;
   - `crash`, `oom`, `timeout`, `undecided`: the bracket is left alone
     and the campaign stops for a human decision.

   Every step is appended to `bisection_history.json` in the campaign
   directory, which also makes the campaign resumable
   (`make bisection-resume`).
2. **Monitoring.** Resource monitoring and ETA estimation are mandatory
   for every run:
   - memory: the runner samples the container's cgroup v2 accounting and
     aborts a run above a fraction of the container limit (90% by
     default). The container has swap disabled, so the kernel OOM killer
     is the backstop;
   - ETA: per run, from the simulation-time rate in the Carpet output
     table; per campaign, from the remaining steps and the wall time of
     the completed runs.

   Both go into `status.json`; `make status` and `make monitor` read it.
3. **Analysis** (`scripts/analyze_critical_exponent.py`). Linear fit of
   `ln M_AH` against `ln(A - A*)` over the supercritical runs. It prints
   `gamma` and `R^2` and writes `critical_scaling.png`. `A*` can be the
   last subcritical bound, the bracket midpoint, an explicit value, or a
   free parameter of a nonlinear fit.

## Directory layout

```text
.
├── Dockerfile                  # overlay on gw230529-et:local (OMP_NUM_THREADS, PATH)
├── docker-compose.yml          # container: CPU set, memory limit without swap, mounts
├── .env.example                # local settings template (copy to .env)
├── Makefile                    # `make help` lists every target
├── makefiles/
│   ├── docker.mk               # image and container lifecycle, docker-check
│   ├── sim.mk                  # single-run, bisection, dry-run
│   ├── monitor.mk              # status, monitor
│   ├── analyze.mk              # analyze
│   └── test.mk                 # check, test
├── templates/
│   └── brill_collapse.par.template   # ET parfile with __AMPLITUDE__ / __FINAL_TIME__
├── par/
│   └── mpi_check.par           # MPI communicator sanity check for docker-check
├── scripts/
│   ├── bisection_runner.py     # bisection driver
│   ├── run_monitor.py          # memory guard, ETA, `status` command
│   ├── cactus_outputs.py       # readers for AHFinderDirect / Carpet output
│   └── analyze_critical_exponent.py  # gamma fit and plot
└── tests/                      # pytest suite (runs inside the container)
```

Campaign output goes to `SIM_OUTPUT_DIR` on the host (default
`~/brill-collapse-output`), mounted at `/home/etuser/simulations`:

```text
<campaign>/
├── bisection_history.json      # every step: A, outcome, M_AH, wall time, peak memory
├── status.json                 # live status for `make status`
└── iter00_A5.0/
    ├── brill_collapse.par      # rendered parfile
    ├── cactus.log              # mpirun stdout and stderr
    ├── resources.csv           # memory samples
    └── output/                 # Cactus output (BH_diagnostics.ah1.gp, h.t*.ah1.gp, ...)
```

## Open issues

- **Initial data.** The stock `IDBrillData` thorn implements exactly the
  Holz q function, but solves the Hamiltonian constraint with the
  PUGH-era `EllSOR` solver, which only works on a unigrid. On a Carpet
  refinement hierarchy it would solve every level independently with its
  own outer boundary condition. A unigrid run at `dx = 0.025` is out of
  reach, so a Carpet-compatible solver (or an imported 2D solution) is
  needed.
- **Gauge near the threshold.** Moving-puncture gauge is reported to
  fail for near-critical Brill data
  ([Hilditch et al. 2013](https://arxiv.org/abs/1309.5008)), so crashes
  close to `A*` are expected.
- **Off-centre horizons.** Near the threshold this family forms disjoint
  apparent horizons around curvature peaks on the symmetry axis, away
  from the origin
  ([Hilditch et al. 2017](https://arxiv.org/abs/1706.01829)). A single
  origin-centred horizon finder can miss them.
- **Universality.** Recent multi-code work finds no universal critical
  solution for vacuum collapse
  ([Baumgarte et al. 2023](https://doi.org/10.1103/PhysRevLett.131.181401)),
  so a measured `gamma` is a property of this family and this setup,
  not a universal constant.

## License

GPL-2.0-or-later; see [LICENSE](LICENSE). The Einstein Toolkit
components used at run time are distributed under their own licenses.
