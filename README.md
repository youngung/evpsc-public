# ΔEVPSC: Elasto-Visco-Plastic Self-Consistent polycrystal model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23999006.svg)](https://doi.org/10.5281/zenodo.23999006)

ΔEVPSC (delta-EVPSC) is a Fortran implementation of an
elasto-visco-plastic self-consistent (EVPSC) polycrystal plasticity
model based on a single homogeneous effective medium (HEM). It can be run as a stand-alone
program or used as a material subroutine in the Abaqus finite element
solver.

EVPSC is derived from the visco-plastic self-consistent (VPSC) code of
R. A. Lebensohn and C. N. Tome (Los Alamos National Laboratory). See
[NOTICE](NOTICE) for details.

## Features

* Stand-alone EVPSC simulations under mixed stress/strain boundary
  conditions [1, 2]
  - A single HEM accounts for both elasticity and visco-plasticity,
    which is computationally faster than EVPSC formulations that use
    two separate HEMs [1, 3]
  - Slip and twinning; Voce and dislocation-density-based hardening
  - Texture and grain-shape evolution
* Abaqus/Standard UMAT interface (`umat/src/evpsc_incr`) that shares
  the core subroutines with the stand-alone program [2]
  - **Note:** the Abaqus/Explicit VUMAT (`vumat.f`) is not complete and
    does not work yet. It is included only as a placeholder for future
    development.
* Forming limit diagram (FLD) prediction based on the
  Marciniak-Kuczynski approach [4]

## Repository layout

| Path                   | Contents                                         |
|:-----------------------|:-------------------------------------------------|
| `src/vpsc_src/`        | Fortran source of the stand-alone EVPSC program and the single-crystal file parser (`parser_fnsx.f`) |
| `umat/src/evpsc_incr/` | Abaqus UMAT interface (VUMAT: placeholder only)  |
| `umat/src/pheno/`      | Phenomenological constitutive models for Abaqus  |
| `pyscripts/`           | Helper scripts (`gen_dim.py`, `build.py`, ...)   |
| `*.mk`, `Makefile`     | Build rules                                      |

## Requirements

* A Fortran compiler: GNU Fortran (`gfortran`) is used for development
  of the stand-alone program; the Intel Fortran compiler is used for
  Abaqus simulations.
* GNU Make
* Python 3 with `numpy`, `matplotlib` and `chardet`. The `python`
  command must start Python 3, because the scripts call each other
  through `python`.

  ```sh
  pip install numpy matplotlib chardet
  ```

* For plotting and some options of the Abaqus driver
  (`umat/src/evpsc_incr/run.py`), the Python packages `MP`
  ([mpl-lib3](https://github.com/youngung/mpl-lib3)) and `TX`
  ([texture3](https://github.com/youngung/texture3))
* Abaqus/Standard for finite element simulations

## Build the stand-alone program

1. Generate the dimension file `src/vpsc_src/parameters.dim`, which
   sets the maximum array sizes (phases, grains, deformation modes,
   etc.):

   ```sh
   python pyscripts/gen_dim.py -o src/vpsc_src/parameters.dim \
       --nphmx 1 --ngrmx 2000 --nmodmx 2 --ntwmmx 0 --nsysmx 24 --ntwsmx 0
   ```

   Run `python pyscripts/gen_dim.py -h` for all options. Alternatively,
   copy one of the provided examples:

   ```sh
   cp src/vpsc_src/parameters_rex.dim src/vpsc_src/parameters.dim
   ```

2. Compile:

   ```sh
   make -j evpsc
   ```

   This creates the executable `evpsc`. Use `make evpsc ENV=debug` for a
   build with run-time checks.

The input files follow the VPSC conventions (`vpsc7.in`, single-crystal
and texture files). For a description of the VPSC input format, see
Chapter 5 of Tome and Lebensohn (2023) [6].

## Run the Abaqus example

A one-element Abaqus/Standard example is in `umat/src/evpsc_incr/`.
After generating `src/vpsc_src/parameters.dim` (step 1 above):

```sh
cd umat/src/evpsc_incr
python run.py --fn umat_example.in --ngrs 100 --cpus 1
```

See [umat/src/evpsc_incr/README.md](umat/src/evpsc_incr/README.md) for
the input file format, the example files and the options of `run.py`.

## License

EVPSC is released under the BSD 3-Clause License; see [LICENSE](LICENSE).
Portions derived from VPSC are (c) 2022 Triad National Security, LLC,
also under the BSD 3-Clause License; see [NOTICE](NOTICE).

## How to cite

If you use EVPSC, please cite Ref. [1] and, for the finite element
interface, Ref. [2]. Citation metadata is in
[CITATION.cff](CITATION.cff).

## References

1. Y. Jeong, C. N. Tome, An efficient elasto-visco-plastic
   self-consistent formulation: Application to steel subjected to
   loading path changes, International Journal of Plasticity 135 (2020)
   102812. https://doi.org/10.1016/j.ijplas.2020.102812
2. Y. Jeong, B. Jeon, C. N. Tome, Finite element analysis using an
   incremental elasto-visco-plastic self-consistent polycrystal model:
   FE simulations on Zr and low-carbon steel subjected to bending,
   stress-relaxation, and unloading, International Journal of
   Plasticity 147 (2021) 103110.
   https://doi.org/10.1016/j.ijplas.2021.103110
3. B. Jeon, Y. Jeong, C. N. Tome, A critical discussion of
   elasto-visco-plastic self-consistent (EVPSC) models, Journal of
   Materials Research and Technology 33 (2024) 7596-7609.
   https://doi.org/10.1016/j.jmrt.2024.11.043
4. Y. Jeong, M.-S. Pham, M. Iadicola, A. Creuziger, T. Foecke, Forming
   limit predictions using a self-consistent crystal plasticity
   framework: a case study for body-centered cubic materials, Modelling
   and Simulation in Materials Science and Engineering 24 (2016) 055005.
   https://doi.org/10.1088/0965-0393/24/5/055005
5. R. A. Lebensohn, C. N. Tome, A self-consistent anisotropic approach
   for the simulation of plastic deformation and texture development of
   polycrystals: Application to zirconium alloys, Acta Metallurgica et
   Materialia 41 (1993) 2611-2624.
   https://doi.org/10.1016/0956-7151(93)90130-K
6. C. N. Tome, R. A. Lebensohn, Material Modeling with the
   Visco-Plastic Self-Consistent (VPSC) Approach: Theory and Practical
   Applications, Elsevier (2023).

## Contact

Youngung Jeong, Department of Materials Science and Metallurgy,
Changwon National University, Republic of Korea
([youngung.github.io](https://youngung.github.io))
