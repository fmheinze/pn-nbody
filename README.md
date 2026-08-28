# pn-nbody
Software for the direct integration of the Newtonian and post-Newtonian equations of motion of the gravitational N-body problem (for more information, see [Heinze, Schäfer & Brügmann 2026](https://arxiv.org/abs/2602.06961)).

## Features
* Direct integration of nonspinning post-Newtonian N-body dynamics through conservative 2PN order, including the leading dissipative 2.5PN radiation-reaction terms
* Complete general-N 2PN dynamics, including the numerical evaluation of the four-body $U_{\mathrm{TT}}^{(4)}$ contribution for $N \geq 4$
* Multiple integration methods: the adaptive Cash–Karp embedded Runge–Kutta RK5(4) method, the second-order symmetric and symplectic implicit midpoint method, and the classical fourth-order Runge–Kutta method
* Configurable initial-condition presets for binaries, hierarchical triples, binary–single and binary–binary scattering, figure-eight configurations, Newtonian Plummer clusters, and relativistic monoenergetic clusters
* Automatic merger handling with active-body bookkeeping, remnant-generation tracking, and merger-event output
* NR-calibrated analytic prescriptions for the remnant mass, spin, and recoil kick, including Lousto–Zlochower and Barausse-based options
* Optional impulse splitting to reduce the number of expensive $U_{\mathrm{TT}}^{(4)}$ evaluations
* Parallel evaluation with OpenMP when available: over four-body subsets for the numerical $U_{\mathrm{TT}}^{(4)}$ contribution, and over the $O(N^3)$ three-body 2PN terms, which dominate the equations of motion whenever $U_{\mathrm{TT}}^{(4)}$ is not in use. Both are arranged so that results do not depend on the number of threads
* Python tools for static trajectory plots and two- or three-dimensional animations

To our knowledge, `pn-nbody` is the first publicly available code to implement the complete general-$N$ dynamics through 2PN order, including the four-body $U_{\mathrm{TT}}^{(4)}$ contribution.

## Todo
* Further performance optimizations (e.g. approximations for objects further away)
* Regularization / handle close encounters
* Additional ODE integration methods (e.g. higher-order symplectic methods, Hermite integrators)
* Spins and spin evolution

## Requirements
* C compiler (e.g. GCC or Clang)
* Linux or macOS
* Optional: OpenMP for parallel $U_{\mathrm{TT}}^{(4)}$ evaluation
* Optional: Python with NumPy and Matplotlib (for visualizations)

## Building
Clone the repository and run `make`:

```bash
git clone https://github.com/fmheinze/pn-nbody.git
cd pn-nbody
make
```

The Makefile automatically checks whether the selected compiler has a working OpenMP setup. If OpenMP is found, it is enabled automatically; otherwise the code is built serially. On macOS with Apple Clang, Homebrew `libomp` is detected automatically when installed:

```bash
brew install libomp
```

The OpenMP behavior can be controlled explicitly:

```bash
make USE_OPENMP=auto   # default: use OpenMP if a working setup is found
make USE_OPENMP=1      # require OpenMP; fail if it is unavailable
make USE_OPENMP=0      # force a serial build
```

A different compiler can be selected through `CC`, for example:

```bash
make CC=gcc
make CC=gcc-15
```

Custom compiler flags can also be supplied:

```bash
make CFLAGS="-O2 -g"
```

## Quick Start Guide
1. After building the executable, you should be able to run the program. You can try one of the examples provided in the `./test` directory, e.g.:
```
./exe/pn-nbody ./test/test_figure_eight.par
```

2. After the simulation has finished, the output should be located in the `./output` directory. You can visualize the output with the included Python scripts in the `./viz` directory, e.g.:
```
python3 ./viz/trajectory_anim.py ./output/test_figure_eight/output_pos.dat
```

3. You can check out the other files in the `./test` directory to see more examples of what is possible.

### Selecting output

Use the `output` parameter to select which output files are created. Available quantities are
`mass`, `position`, `momentum`, `velocity`, `spin`, `energy`, and `merger`; separate them with
spaces:

```text
output = position velocity energy merger
```

This example creates only `output_pos.dat`, `output_vel.dat`, `output_energy.dat`, and
`output_merger.dat`. If the parameter is omitted, its default selects every available quantity.

To monitor a selected pair, add `orbit_X_Y` to `output`, where `X` and `Y` are distinct 1-based
initial body numbers. A parameter with the same name selects the columns in that pair's file:

```text
output = position energy orbit_1_2 orbit_1_3
orbit_1_2 = semimajor_axis_newtonian eccentricity_newtonian semimajor_axis_radial eccentricity_radial
orbit_1_3 = eccentricity_newtonian
```

This creates `output_orbit_1_2.dat` and `output_orbit_1_3.dat`. The Newtonian quantities are
instantaneous osculating elements computed from the relative position and the full coordinate
velocity at all enabled PN orders. `semimajor_axis_newtonian` uses the pair's Newtonian specific
orbital energy (and is negative for a hyperbolic orbit); `eccentricity_newtonian` is the magnitude
of its Newtonian eccentricity vector.

The relativistic radial quantities are `pericenter_radial`, `apocenter_radial`,
`semimajor_axis_radial`, and `eccentricity_radial`. For each selected pair, the code removes its
center-of-mass momentum with the canonical Jacobi decomposition and evaluates the isolated-pair
ADM Hamiltonian using the enabled conservative 1PN and 2PN terms (on top of the Newtonian term).
It then solves `H(r, p_r=0, J) = E` for the two turning points enclosing the current separation and
defines

```text
semimajor_axis_radial = (pericenter_radial + apocenter_radial) / 2
eccentricity_radial   = (apocenter_radial - pericenter_radial)
                        / (apocenter_radial + pericenter_radial)
```

Other bodies and the dissipative 2.5PN term are excluded from this instantaneous, osculating
definition. Unbound states and states for which two positive enclosing turning points cannot be
found produce `nan` for all radial quantities. The four radial columns share a single root solve;
requesting only radial quantities does not compute coordinate velocities.

Each endpoint follows its merger lineage: if it merges with a third body, subsequent rows use the
remnant. If the two monitored lineages merge with each other, subsequent orbital values are `nan`.

## License

Copyright (c) 2026, Felix M. Heinze

This project is licensed under the MIT License. See the LICENSE file for details.
