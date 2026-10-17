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

> **Status: Phase 1b done, campaign tooling in place.** Runs are
> monitored (memory guard, ETA), and the evolution, mesh refinement and
> horizon finder work in axisymmetry with `Cartoon2D` (Phase 1a). The
> Brill initial data are solved in 2D and loaded onto the grid by a small
> thorn; the ADM masses reproduce every published value for this family.
> The three-way classifier, the resumable bisection and the two-sided
> gamma fit are implemented and tested on synthetic output; their
> thresholds await calibration. The gauge test at `A = 5` (Phase 2) is
> next.

## Plan

| Phase | Content |
| --- | --- |
| 1a | Monitored single runs (memory guard, ETA). `Cartoon2D` feasibility with Brill-Lindquist black holes: ML_BSSN + mesh refinement + time refinement + AHFinderDirect, compared with the 3D octant grid |
| 1b | Initial data: a Python solver for the axisymmetric Hamiltonian constraint, interpolated onto every Carpet level by a small thorn; checked against the published ADM mass (0.70 at `A = 5`) |
| 2 | Gauge feasibility at `A = 5`: reproduce the reported 1+log failure, then try harmonic and shock-avoiding slicing |
| 2.5 | Coarse pilot campaign: horizon positions, classifier thresholds, rough `A*` |
| 3 | Freeze the numerical setup; three-way classifier (collapse / dispersion / undecided); known-answer tests |
| 4 | Production bisection |
| 5 | Scaling fits on both sides: `M_AH` above `A*`, maximum curvature below |

Success criterion: with `Cartoon2D`, type II behaviour near the
threshold (`M_AH → 0` above `A*`, growing peak curvature below it);
without it, a robust bracket on `A*`. A precise `gamma` is optional.

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
  which this project overlays. If it exists under another tag (for
  example a registry tag), `docker tag` it to `gw230529-et:local` or set
  `BASE_IMAGE` in `.env`. Otherwise build it once in
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

Initial data: solve the constraint for one amplitude and write the table
that the `BrillPsiTable` thorn loads (`SIM_OUTPUT_DIR/tables/`):

```bash
make id-table A=5.0
make run-par PAR=par/initial_data/brill_id_check.par CAMPAIGN=initial_data
```

After editing a thorn under `thorns/`, `make thorn-build` recompiles it
in the running container; `make docker-build` bakes it into the image.

Any parfile can run under the monitor, for example a feasibility test:

```bash
make run-par PAR=par/feasibility/cartoon_bl_single.par CAMPAIGN=feasibility
make status CAMPAIGN=feasibility
```

Campaign (once the production template is in place, Phase 3):

```bash
make single-run AMP=4.0                         # one classified run
make bisection EXTRA="--verify-endpoints"       # check both ends, then bisect (A_MIN, A_MAX, MAX_ITER, CAMPAIGN)
make status                                     # current amplitude, run and campaign ETA, memory
make bisection-resume EXTRA="--set-outcome 3=collapse"   # after a stop: record a human verdict
make analyze                                    # gamma on both sides, reports/critical_scaling.png
make reclassify RUN=/home/etuser/simulations/<campaign>/single/iter00_A4.8   # after a classifier change
```

Figures for the wiki (rendered into `reports/figures/`, which is gitignored):

```bash
make figures                                    # initial data and Phase 2 gauge figures
make wiki-figures                               # copy them into ../brill-critical-collapse-einstein-toolkit.wiki/images/
```

## Workflow

1. **Bisection** (`scripts/bisection_runner.py`). For each step the
   runner renders the template at the bracket midpoint (solving the
   initial data if the template asks for a table), runs `cactus_sim`
   under `mpirun.mpich`, and classifies the run
   (`scripts/outcome_classifier.py`):
   - `collapse`: a horizon found at consecutive find times with a mass
     that does not shrink; the upper bound moves down;
   - `dispersion`: positive evidence only: no horizon ever, the peak
     curvature invariant far below its initial value for a while and
     not growing back, and a recovered lapse; the lower bound moves up;
   - `crash`, `oom`, `timeout`, `undecided`: the bracket is left alone
     and the campaign stops for a human decision. Missing a horizon and
     calling dispersion too early would both bias `A*` upwards, so the
     classifier abstains rather than guesses.

   Runs stop early once decided (horizon mass settled, or dispersion
   established). Every run is appended to `bisection_history.json` with
   its reasons and evidence, together with a fingerprint of the
   numerical setup; resuming with a different template, classifier
   thresholds or executable is refused. `--set-outcome ITER=...`
   records a human verdict for a stopped run, `--retry` reruns it.
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
3. **Analysis** (`scripts/analyze_critical_exponent.py`). Fits
   `ln M_AH` against `ln(A - A*)` over the collapse runs and the peak
   curvature `ln K_max` against `ln(A* - A)` over the dispersion runs
   (`K_max ~ (A* - A)^(-4 gamma)`), so `gamma` is measured twice,
   independently. It prints both with their errors and writes
   `critical_scaling.png`. `A*` can be the last subcritical bound, the
   bracket midpoint, an explicit value, or a free parameter of a
   nonlinear fit.

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
│   ├── figures.mk              # figures, wiki-figures
│   └── test.mk                 # check, test
├── templates/
│   ├── brill_collapse.par.template   # 3D octant parfile with __AMPLITUDE__ / __FINAL_TIME__
│   └── brill_cartoon.par.template    # Cartoon2D parfile used from Phase 2 on
├── par/
│   ├── mpi_check.par           # MPI communicator sanity check for docker-check
│   ├── feasibility/            # Phase 1a Cartoon2D tests with Brill-Lindquist black holes
│   ├── initial_data/           # Phase 1b checks of the Brill initial data on the grid
│   └── gauge/                  # Phase 2 slicing / shift variants of the Cartoon2D template
├── patches/                    # source patches to ET thorns, applied in the overlay build
├── thorns/
│   └── BrillCollapse/BrillPsiTable/  # loads the psi table onto the Cactus grid
├── scripts/
│   ├── bisection_runner.py     # bisection driver, history, resume
│   ├── outcome_classifier.py   # collapse / dispersion / undecided
│   ├── run_monitor.py          # memory guard, ETA, `status` command
│   ├── cactus_outputs.py       # readers for AHFinderDirect / Carpet output
│   ├── brill_initial_data.py   # 2D Hamiltonian constraint solver, psi table writer
│   ├── refine_template.py      # add finer refinement levels to a parfile template
│   ├── apply_cactus_patches.sh # idempotent application of patches/
│   ├── analyze_critical_exponent.py  # gamma fit and plot
│   ├── plot_initial_data.py    # figures of the initial data (ADM mass, psi)
│   ├── plot_runs.py            # figures comparing finished runs (horizon mass, lapse, curvature)
│   └── plot_style.py           # shared matplotlib style of the figures
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
  PUGH-era `EllSOR` solver, which only works on a unigrid. This project
  solves the axisymmetric constraint in Python
  (`scripts/brill_initial_data.py`, 6th-order, compactified polar grid)
  and the `BrillPsiTable` thorn interpolates it onto the grid. The ADM
  masses agree with [Hilditch et al. 2017](https://arxiv.org/abs/1706.01829)
  and [Khirnov & Ledvinka 2019](https://arxiv.org/abs/1908.06034) to all
  quoted digits; the oblate masses quoted in
  [Hilditch et al. 2013](https://arxiv.org/abs/1309.5008) (Sec. IV.3) do
  not match and look like a transcription error.
- **Gauge near the threshold.** 1+log slicing fails for near-critical
  Brill data before a horizon forms: at `A = 5` (`A* ≈ 4.7`) the slices
  become non-smooth at `t ≲ 4.9` at any resolution
  ([Khirnov & Ledvinka 2019](https://arxiv.org/abs/1908.06034); see also
  [Hilditch et al. 2013](https://arxiv.org/abs/1309.5008)). Harmonic and
  shock-avoiding slicing
  ([Baumgarte & Hilditch 2022](https://arxiv.org/abs/2207.06376)) are
  the planned fallbacks.
- **Off-centre horizons.** Near the threshold the global curvature peak
  sits on the symmetry axis, away from the origin
  ([Hilditch et al. 2017](https://arxiv.org/abs/1706.01829);
  [Suárez Fernández et al. 2022](https://arxiv.org/abs/2205.04379)), so
  horizons can form as a pair at `(0, 0, ±z0)`. A single origin-centred
  horizon finder can miss them; on-axis finders are added.
- **Universality.** Recent multi-code work finds no universal critical
  solution for vacuum collapse
  ([Baumgarte et al. 2023](https://doi.org/10.1103/PhysRevLett.131.181401)),
  so a measured `gamma` is a property of this family and this setup,
  not a universal constant.

## License

GPL-2.0-or-later; see [LICENSE](LICENSE). The Einstein Toolkit
components used at run time are distributed under their own licenses.
