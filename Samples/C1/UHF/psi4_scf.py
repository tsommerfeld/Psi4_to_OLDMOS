#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 12:50:17 2022

@author: Thomas Sommerfeld
"""


import os
import sys
curr_dir=os.getcwd()
p=curr_dir.find('Samples')
root=curr_dir[:p]
sys.path.append(root+'lib')

import argparse

import psi4
import c4_comp_geo
import P4toC4_aux

def main():

    parser = argparse.ArgumentParser(description='Psi4 SCF')
    parser.add_argument("fname", type=str, 
        help="Cfour output file for computational symmetry and orientation")
    parser.add_argument("-b",  type=str, default="STO-3G", help='basis set')
    parser.add_argument("-v", type=int, default=1, help='verbosity')
    arguments = parser.parse_args()

    verbose = arguments.v
    Basis = arguments.b         
    sym, coors = c4_comp_geo.cfour_comp_sym_and_geo(arguments.fname)

    if verbose > 0:
        print('Molecule from', arguments.fname)
        print(coors)
        print('Basis is', Basis)
        
    psi4.set_memory('500 MB')
    psi4.core.set_global_option("BASIS", Basis)
    psi4.core.set_global_option("SCF_TYPE", "pk")
    psi4.core.set_global_option("REFERENCE", "UHF")
    psi4.core.set_global_option("D_CONVERGENCE", 1e-8)
    psi4.core.set_global_option("PUREAM", "True")
    psi4.core.set_output_file('output.dat', False)

    mol_str =  '0 2\n' + coors + 'symmetry ' + sym + '\n'
    mol_str += 'units Bohr\nno_reorient'
    if verbose > 1:
        print('-------------------')
        print(mol_str)
        print('-------------------')
    mol = psi4.geometry(mol_str)
    E, wf = psi4.energy('scf', return_wfn=True, molecule=mol)
    print(f'HF energy = {E}')

    P4toC4_aux.make_OLDMOS(wf, verbose=verbose)

    if verbose > 0:
        print('Psi4-MOs written to PSIMOS')
    return

if __name__ == "__main__":
    main()

