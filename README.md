# pn-nbody
Software for the direct integration of the Newtonian and post-Newtonian equations of motion of the gravitational N-body problem (for more information, see [Heinze, Schäfer & Brügmann 2026](https://arxiv.org/abs/2602.06961)).

## Features
* Direct integration of the post-Newtonian equations of motion up to 2.5PN order, suitable for relativistic simulations of few-body systems (first code to fully include all post-Newtonian terms for N bodies at 2PN order)
* Implemented methods: adaptive embedded 5th-order Runge-Kutta method, the symmetric and symplectic implicit midpoint method, and the 4th-order Runge-Kutta method
* Convenient use of customizable initial condition presets (e.g. binary, binary-single scattering, binary-binary scattering, figure-eight orbit)

## Todo
* Performance optimizations (parallelization, approximations for larger N)
* Regularization / handle close encounters
* Implementation of additional ODE integration methods (e.g. higher-order symplectic methods, Hermite integrators)

## Requirements
* C compiler (e.g. gcc/clang)
* Linux/macOS
* Optional: Python with numpy and matplotlib (for visualizations)
* Optional: Cuba 4.2.2 (If you want to solve the full 2PN equations of motion for N bodies, with N > 3, it is necessary to compute the four-point correlation function UTT4, see e.g. [Heinze, Schäfer & Brügmann 2026](https://arxiv.org/abs/2602.06961). UTT4 contains an integral that can currently only be solved numerically, which is computationally very expensive. The numerical evaluation requires the [Cuba library](https://feynarts.de/cuba/). If you want to do that, download and install it. You can compile and run the program also without the Cuba library, but then you cannot include UTT4.)

## How To Install
### Option 1: Build without the Cuba library
```bash
git clone https://github.com/fmheinze/pn-nbody.git
cd pn-nbody
make
```

### Option 2: Build with the Cuba library
```bash
git clone https://github.com/fmheinze/pn-nbody.git
cd pn-nbody
make CUBA_DIR=/path/to/Cuba
```

## Quick Start Guide
1. After building the executable, you should be able to run the program. You can try one of the examples provided in the `./test` directory, e.g.:
```
./exe/pn-nbody ./test/test_figure_eight.par
```

2. After the simulation has finished, the output should be located in the `./output` directory.  You can visualize the output with the included python scripts in the `./viz` directory, e.g.:
```
python3 ./viz/trajectory_anim.py ./output/test_figure_eight/output_pos.dat
```

3. You can check out the other files in the `./test` directory to see more examples of what is possible.

## License

Copyright (c) 2026, Felix M. Heinze

This project is licensed under the MIT License. See the LICENSE file for details.