# EVPSC-FE: Abaqus/Standard UMAT

This folder contains the Abaqus/Standard user material subroutine
(UMAT) of EVPSC and the Python driver `run.py`, which builds the UMAT
and runs an Abaqus job. The model and its implementation are described
in Jeong, Jeon and Tome, International Journal of Plasticity 147 (2021)
103110.

**Note:** the Abaqus/Explicit VUMAT (`vumat.f`) is not complete and does
not work yet. It is included only as a placeholder for future
development.

## Requirements

* Abaqus/Standard with a Fortran compiler supported by your Abaqus
  version (Intel Fortran has been used by the developers)
* Python 3 with `numpy`, `matplotlib` and `chardet`
* For plotting and the `--j2`/`--elastmod` options: `MP`
  ([mpl-lib3](https://github.com/youngung/mpl-lib3)) and `TX`
  ([texture3](https://github.com/youngung/texture3))

## Files needed for a simulation

An EVPSC-FE simulation needs four files:

1. A texture file (discrete grains, or an ODF in MTEX format)
2. A single-crystal file: elastic constants, slip/twin systems and
   hardening parameters. The formats of the texture and single-crystal
   files are the same as in VPSC7/VPSC8.
3. An Abaqus input file (`*.inp`) that uses the user material
4. An EVPSC-FE input file that lists the three files above together
   with a few flags

## Example

This folder includes a one-element example:

| File                    | Contents                                         |
|:------------------------|:-------------------------------------------------|
| `umat_example.in`       | EVPSC-FE input file                              |
| `UMAT_DEV_VEL_0001.inp` | Abaqus input: one C3D8 element, tension then compression |
| `ferrite_odf.odf`       | Ferrite ODF (MTEX format)                        |
| `F.sx`                  | Ferrite single-crystal file                      |
| `bcc.dif`               | Diffraction planes (used only if `idiff=1`)      |
| `texinlc.tex`           | Alternative: low-carbon steel texture (discrete) |
| `lc_fstar7.sx`, `lc.txt`| Alternative: low-carbon steel single-crystal files |

To run it:

1. Generate the dimension file `src/vpsc_src/parameters.dim` (see the
   top-level [README](../../../README.md)).
2. From this folder, run:

   ```sh
   python run.py --fn umat_example.in --ngrs 100 --cpus 1
   ```

   `--ngrs` is the number of grains per phase at each integration
   point. Run `python run.py --help` for all options.

`run.py` creates a new working folder, builds the UMAT there and starts
the Abaqus job.

**Note:** on Linux, `run.py` copies the `abaqus_v6.env` file from the
installation path of the cluster used by the developers (see `run.py`).
Adjust that path, and the `ABQ` command at the top of `run.py`, for
your system.

## EVPSC-FE input file

`umat_example.in` has the following format. Lines starting with `*` are
labels; the line below each label holds its value.

```text
** input to run dEVPSC-FE run via UMAT
* isolver (0, 1)
0
* dimension file
../../../src/vpsc_src/parameters.dim
* nph
1
* wph
1.0
* iodf (0: discrete texture; 1: odf (mtex format))
1
* texture/ODF file
ferrite_odf.odf
* sx file
F.sx
* idiff (0: no diffraction; 1: diffraction)
0
* diff file (dummy if idiff = 0 is given in the above line)
bcc.dif
* abaqus input file
UMAT_DEV_VEL_0001.inp
* noel and npt subjected to 'lookup' (could be dummy unless <isa> flag is on)
1
1 1
* ipx
0
* pxfn (dummy if ipx is 0)
dummy
* implicit_update (0: False, 1: True)
0
```

To use the low-carbon steel data instead, set `iodf` to `0`, the
texture file to `texinlc.tex` and the sx file to `lc_fstar7.sx`.

## Main options of `run.py`

| Option                   | Description                                                  |
|:-------------------------|:-------------------------------------------------------------|
| `--fn FN`                | EVPSC-FE input file                                          |
| `--cpus CPUS`            | Number of CPUs for Abaqus                                    |
| `--ngrs NGRS [NGRS ...]` | Number of grains per phase at each integration point         |
| `--inp INP`              | Abaqus input file (overrides the one in the input file)      |
| `--iregime IREGIME`      | 1: VPSC; 2: delta-EVPSC                                      |
| `--idd`                  | Use the dislocation-density (DD) hardening model             |
| `--irgvb`                | Use the RGVB dislocation-density hardening model             |
| `--isa`                  | After the FE run, re-run the stand-alone code at a material point |
| `--pathwork PATHWORK`    | Working folder for the FE run                                |
| `--usermat USERMAT`      | Your own user material subroutine instead of `umat.f`        |
| `--j2`                   | J2 plasticity with a flow curve from stand-alone EVPSC       |
| `--elastmod`             | Calculate the elastic moduli and exit                        |
