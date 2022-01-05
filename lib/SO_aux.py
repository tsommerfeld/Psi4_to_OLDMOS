#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 12:13:12 2022

@author: Thomas Sommerfeld

Class for dealing with Psi4 to Cfour SO mapping/transformation

"""

import numpy as np

class SymOrbs:
    """ 
        SOs for one irrep 
        each SO is a column; each column has the length of all AOs
        this matrix is needed to transform one block of MOs
        from the AO to the SO basis
        
        But each column contains at most 8 entries (max depends on group)
        
        To reorder SOs Psi4 -> Cfour, store by index, coefficient 
        
    """

    def __init__(self, C, order=8):
        """ Instance from one symmetry-block of SOs:  nAOs*nSOs[irrep] """
        n, m = C.shape
        self.nbfs = n     # number of AOs in the basis set
        self.nsos = m     # number of SOs in this irrep or block
        self.naos = np.zeros(m, int) # AOs contributing to this SO 
        self.jao  = np.zeros((order,m), int) # indices of the contributing AOs
        self.cao  = np.zeros((order,m)) # coefficients of the contributing AOs

        for jso in range(self.nsos):
            so_vec = C[:,jso]
            nonz = np.nonzero(so_vec)[0]
            n = len(nonz)
            self.naos[jso] = n
            self.jao[:n,jso] = nonz
            self.cao[:n,jso] = so_vec[self.jao[:n,jso]]

    def print(self):
        """ print SOs as linear combinations """
        for jso in range(self.nsos):
            print(f'{jso:2d}', end='')
            for i in range(self.naos[jso]):
                i_ao = self.jao[i,jso]
                c_ao = self.cao[i,jso]
                print(f'  {c_ao:8.5f}({i_ao:3d})', end='')
            print()

    def first_AOs(self):
        """ returns the indices of the first AOs contributing to each SO """
        return self.jao[0,:]

    def inv_coef(self):
        """ 
        return the inverse of the SO vector elements 
        this factor is needed for Cfour-MO = Psi4-MO/inv_coef
        performed by P4toC4_aux.psi4_to_c4()
        """ 
        return 1.0/abs(self.cao[0,:])

    def matrix(self):
        """ recreate matrix representation C[iAO,jSO] in Psi4 order """
        nao, nso = self.nbfs, self.nsos        
        C = np.zeros((nao,nso))
        for jso in range(nso):
            for k in range(self.naos[jso]):
                iao = self.jao[k,jso]
                C[iao,jso] = self.cao[k,jso]
        return C
    
    def cfour_matrix(self, so_p2c):
        """ 
        matrix representation C[iAO,jSO] in Cfour order 
        so_p2c is the psi4_to_cfour SO mapping
        so Cfour_SO[:,i] = Psi4_SO[:,so_p2c[i]]
        
        This is untested in probably not needed.
        Rather reorder and scale the Psi4-MOs in the SO basis 
        using so_p2c.
        
        The scaling coefficients though. 
        
        """
        nao, nso = self.nbfs, self.nsos        
        C = np.zeros((nao,nso))
        for jso in range(nso):
            for k in range(self.naos[jso]):
                iao = self.jao[k,jso]
                C[iao,jso] = self.cao[k,so_p2c[jso]]
        return C
    
    
    
    
    
    
    
    