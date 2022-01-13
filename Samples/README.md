# Tests for various symmetry/reference combinations.

Find whatever matches your case best and work from there.

## Warnings

**For the scripts to work, the following settings are required:**

1. The Psi4 and Cfour executables must be in your `PATH`.
2. The environment variables `CFOUR_SCRATCH` and  `PSI_SCRATCH` must be set and point to **scratch directories**.

## Workflow

1. Use the `cfour.scf` provided in each directory to run STO-3G with your molecule. The output defines the computational orientation and point groups, and it may be useful for the `OCCUPATION`.
2. Run `psi4_scf.py` provided in this directory with the desired basis set. This "should" run fairly automatically and create the file `PSIMOS`.
3. Re-run `cfour-scf` with the desired basis set. This script creates all the files Cfour needs to read the Psi4-MOs as guess.

**Remark:** It is frequently necessary to set the `OCCUPATION` in the second `cfour-scf` run. Otherwise, Cfour will guess the occupation based on its core-guess even though guess orbitals are supplied, and if the guessed occupation does not match the guessed orbitals, the while guess does not work.  


All three steps can be run in one go using the script `run_me`.

## Starting new projects

At the moment, a lot of user input and attention is required:
* The Psi4 basis sets should be put into the P4.BASIS file
  using the psi4.core basisset.genbas() function.
  This is not strictly required, but as Cfour tends to use general contractions, while Psi4 uses segmented, the MOs are slightly different,
  if native basis sets are employed.
* In both cfour_scf and psi4_scf.py, charge, multiplicity, reference, and possibly occupation must be adjusted to your molecule or ion.
