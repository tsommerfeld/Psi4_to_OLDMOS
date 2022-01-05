# Psi4_to_OLDMOS

Exports Psi4 MOs in CFOUR OLDMOS format 

CFOUR is an excellent many-body package, but sometimes it can be tricky to converge
the SCF. Personally I had trouble with
* non-valence anions (one-center with a (7s7p7d) set of diffuse functions).
* temporary anion calculations (several centers with (2s4p2d) sets of diffuse functions with tight exponents).
* cluster anions such as O2(-)..benzene.

In some of these cases, neither varients of the normal tricks [SCF_MAXCYC, SCF_DAMP, SCF_EXTRAP] nor QRHF_GUESS variations nor QCSCF helped. 
I had some discussions with the developers, but make nor mistake, these systems can be hard to converge in any code.

The Psi4_to_OLDMOS interface adds yet another option: Converge the SCF of your choice in Psi4- provided you can- and use the Psi4 orbitals as a guess in CFOUR.

If everything is done exactly the same in both codes, this works as well as it should: CFOUR converges in one iteration.

Currently, the limitations are:
* RHF and UHF
* C1 symmetry
* up to f-functions

## Workflow

1. Obtain a common orientation
  * Cartesian coordinates that will not reorient in Cfour and Psi4.
  * Or read the computational geometry from Cfour into Psi4 and use no_reorient (see Samples).
2. Run the desired SCF in Psi4.
3. For optimal results, use the Psi4 wavefunction object to create GENBAS (see notebooks). Psi4 uses segmented contraction, CFOUR doesn't. 
4. Create OLDMOS as shown in the provided Jupyter notebooks and Python scripts. 
4. Start CFOUR with ZMAT, GENBAS, OLDMOS, and an empty file JFSGUESS present.

## Repository

* `Samples` : complete comparison of the build-in core guess with the Psi4-guess.
* `CFOUR` : Psi4 basis sets (P4.BASIS) and some experimental scripts with symmetry.  
* `lib` : Python auxiliary functions.
* `dev` : uncommented notebooks used for development.
