#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 12:50:17 2022

@author: Thomas Sommerfeld
"""

import numpy as np
import sys
import argparse

import psi4
sys.path.append('../../../lib')
from P4toC4_aux import basis_mapping, psi4_to_c4, write_oldmos
from c4_comp_geo import cfour_comp_sym_and_geo

def main():

    parser = argparse.ArgumentParser(description='Psi4 SCF')
    parser.add_argument("fname", type=str, 
        help="Cfour output file for computational symmetry and orientation")
    parser.add_argument("-b",  type=str, default="STO-3G", help='basis set')
    parser.add_argument("-v", type=int, default=1, help='verbosity')
    arguments = parser.parse_args()

    verbose = arguments.v
    Basis = arguments.b         
    sym, coors = cfour_comp_sym_and_geo(arguments.fname)

    if verbose > 0:
        print('Molecule from', arguments.fname)
        print(coors)
        print('Basis is', Basis)
        
    psi4.set_memory('500 MB')
    psi4.core.set_global_option("BASIS", Basis)
    psi4.core.set_global_option("SCF_TYPE", "pk")
    psi4.core.set_global_option("REFERENCE", "RHF")
    psi4.core.set_global_option("D_CONVERGENCE", 1e-8)
    psi4.core.set_global_option("PUREAM", "True")
    psi4.core.set_output_file('output.dat', False)

    mol_str =  '0 1\n' + coors + 'symmetry ' + sym + '\n'
    mol_str += 'units Bohr\nno_reorient'
    if verbose > 1:
        print('-------------------')
        print(mol_str)
        print('-------------------')
    h2o = psi4.geometry(mol_str)
    E, wf = psi4.energy('scf', return_wfn=True)
    print(f'HF energy = {E}')

    p2c_map, scale = basis_mapping(wf.basisset(), verbose=0)
    Ca=np.array(wf.Ca())
    Cb=np.array(wf.Cb())
    Ca_C4 = psi4_to_c4(Ca, p2c_map, scale)
    Cb_C4 = psi4_to_c4(Cb, p2c_map, scale)
    write_oldmos('PSIMOS', Ca_C4, Cbs=Cb_C4)
    if verbose > 0:
        print('Psi4-MOs written to PSIMOS')
    return

if __name__ == "__main__":
    main()

