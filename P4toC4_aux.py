#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 17:27:58 2021

@author: thomas

  

"""

import numpy as np
import sys
#import psi4


def psi4_to_c4(Cp4, offset):
    """
    reorder the Psi4_MOs into Cfour order based on offset[]
    Cfour_MO[i+offset[i]] = Psi4_MO[i]

    Parameters
    ----------
    Cp4 : np.array(nbf, nbf) Psi4MOs
    
    offset : np.array(nbf) offset vector

    Returns
    -------
    Cc4:  np.array(nbf, nbf) Cfour4MOs

    """
    nbf = len(offset)
    n, m = Cp4.shape
    if n != nbf or m != m:
        msg = 'Error in psi4_to_c4: inconsistent dimensions %d %d %d' % (nbf, n, m)
        sys.exit(msg)
    
    Cc4=np.zeros((nbf,nbf))
    for i in range(nbf):
        Cc4[i+offset[i],:] = Cp4[i,:]
    return Cc4


def ao_offset(basisset, verbose=0):
    """
    offset vector for Psi4 to Cfour MO mapping
    
    Cfour_MO[i+offset[i]] = Psi4_MO[i]

    for s-fns offset is 0
    for p-shells offset is [+2,-1,-1] : (z,x,y) -> (x,y,z)

    Parameters
    ----------
    basisset : psi4.core.basisset
        Psi4 class for basis sets
        we need .nshells and .shell.am

    Returns
    -------
    the offset vector for the basis set

    """
    if verbose > 0:
        print('Creating offset vector')

    ll = {'s': (1, [0]),
          'p': (3, [2, -1, -1])}
    nbf = basisset.nbf()
    offset = np.zeros(nbf, int)

    ibf = 0
    for ishell in range(basisset.nshell()):
        shell = basisset.shell(ishell)
        am = shell.amchar
        l = shell.am
        if verbose > 0:
            print(f"Shell {ishell}  {am} {l}")
        add, p2c = ll[am]
        offset[ibf:ibf+add] = p2c
        ibf = ibf + add
    return offset
    

def read_oldmos(fname, verbose = 1):
    """
    read a Cfour OLDMOS file
    ASSUMES nAOs = nMOs

    Parameters
    ----------
    fname : input file (str)

    Returns
    -------
    nbas : int number_of_basis_fns
    cs : np.array orital_coeffcients[ao,mo]

    """
    if verbose > 0:
        print('reading orbitals from '+fname)
    file = open(fname)
    lines = file.readlines()
    file.close()
    
    nl = len(lines)
    words=lines[-1].split()
    n_orb_last = len(words)
    if verbose > 0:
        print(f'  {nl} lines, {n_orb_last} orbitals at the bottom')
    if n_orb_last == 4:
        size = nl*4
        nbas = int(np.sqrt(size))
        if nbas*nbas != size:
            print(f'{size} != {nbas}**2')
            sys.exit('TROUBLE: basis size cannot be determined.')
    else:
        i = 1
        while len(lines[-i].split()) < 4:
            i += 1
        nbas = i - 1
    print(f'  nAOs = nMOs = {nbas}')
    
    cs = np.zeros((nbas,nbas))
    l = 0
    for jmo in range(0, nbas, 4):
        for ao in range(0, nbas):
            words = lines[l].split()
            for i in range(len(words)):
                cs[ao,jmo+i] = float(words[i])
            l += 1
    return nbas, cs



def write_oldmos(fname, Cs):
    """
    Write MOs in Cfour OLDMOS format
    That means packages of four MOs   
    C[0,0] C[0,1] C[0,2] C[0,3]
    C[1,0] C[1,1] C[1,2] C[1,3]
    C[2,0] C[2,1] C[2,2] C[2,3]
    ...    ...     ...    ...

    Format for each individual coefficient: 30.20E

    Parameters
    ----------
    fname : str filename
    Cs : np.array with MO coefficients Cs[ao,mo]

    Returns
    -------
    None.
    """

    f = open(fname, 'w')
    nbf = Cs.shape[0]
    for j in range(0, nbf, 4):
        """ Normally we write groups of 4 MOs. The last group may be smaller. """
        ngr = 4
        if j + ngr >= nbf:
            ngr = nbf - j
            
        for ao in range(nbf):
            line = ''
            for igr in range(ngr):
                line += f"{Cs[ao,j+igr]:30.20E}"
            line += "\n"
            f.write(line)
    f.close()
    return