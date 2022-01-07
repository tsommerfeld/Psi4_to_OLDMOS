#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 2022

@author: Thomas Sommerfeld
"""

import os
import sys
curr_dir=os.getcwd()
p=curr_dir.find('Samples')
root=curr_dir[:p]
sys.path.append(root+'lib')

import numpy as np
import argparse

import psi4
from P4toC4_aux import basis_mapping, psi4_to_c4, write_oldmos
from c4_comp_geo import cfour_comp_sym_and_geo
from SO_aux import SymOrbs

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
    mol = psi4.geometry(mol_str)
    E, wf = psi4.energy('scf', return_wfn=True, molecule=mol)
    print(f'HF energy = {E}')

    p2c_map, p2c_scale = basis_mapping(wf.basisset(), verbose=0)
    naos=len(p2c_map)
    c2p_map=np.zeros(naos, int)
    for i in range(naos):
        c2p_map[p2c_map[i]] = i 

    # MOs and SOs
    Ls=wf.aotoso()
    C_SO=wf.Ca()
    if verbose > 0:
        print('SO dimensions:', Ls.shape)
        print('MO dimensions:', C_SO.shape)        

    irrep_lst = []
    
    for isym in range(wf.nirrep()):
        SOs=SymOrbs(Ls.nph[isym], order=wf.nirrep())
        if verbose > 3:
            SOs.print()
        p4_first_AOs = SOs.first_AOs()
        cfour_first_AOs = p2c_map[SOs.first_AOs()]
        ao_scale = p2c_scale[SOs.first_AOs()]
        so_c2p = np.argsort(cfour_first_AOs)
        nsos=len(so_c2p)
        so_p2c=np.zeros(nsos, int)
        for i in range(nsos):
            so_p2c[so_c2p[i]] = i
        so_scale=SOs.inv_coef()
        scale = so_scale*ao_scale
        if verbose > 1:
            print(f'\nIrrep {isym}')
            print('AO-order  AO-order   Cfour    argsort    AO     SO')
            print('  Psi4     Cfour    argsort   inverted  scale  scale')
            for i in range(nsos):
                print(f'{p4_first_AOs[i]:4d}{cfour_first_AOs[i]:9d}', end='')
                print(f'{so_c2p[i]:11d}{so_p2c[i]:10d}', end='')
                print(f'{ao_scale[i]:11.3f}{so_scale[i]:7.3f}')
        
        C=psi4_to_c4(C_SO.nph[isym], so_p2c, scale, use_scale=True)
        irrep_lst.append(C)
                
    C_SOr = psi4.core.Matrix.from_array(irrep_lst)
    for irrep in range(wf.nirrep()):
        mode = 'w'
        if irrep > 0:
            mode = 'a'
        write_oldmos('PSIMOS', C_SOr.nph[irrep], mode=mode)

    if verbose > 0:
        print('Psi4-MOs written to PSIMOS')
    return

if __name__ == "__main__":
    main()

