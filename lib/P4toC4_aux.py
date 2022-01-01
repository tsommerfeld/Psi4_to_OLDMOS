#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 17:27:58 2021

@author: Thomas Sommerfeld

"""

import numpy as np
import sys
import psi4



def AO_MO_by_irrep(wfn, verbose=0):
    """
    creates the psi4.core.Matrix C_AO_pi; each entry in C_AO_pi.nph has
        index 0: AOs
        index 1: MO[irrep]

    Parameters
    ----------
    wfn : psi4.core.wavefunction (so far RHF)

    Returns
    -------
    psi4.core.Matrix C_AO_pi
    
    """
    eps_pi = wfn.epsilon_a()
    eps_by_sym = np.concatenate(eps_pi.nph)
    j_sym2E = np.argsort(eps_by_sym)
    n = len(j_sym2E)
    j_E2sym = np.zeros(n, int)
    for i in range(n):
        j_E2sym[j_sym2E[i]] = i
    if verbose > 1:
        print(' MO  sym-to-E   E-to-sym')    
        for i in range(n):    
            print(f'{i:3d}   {j_sym2E[i]:7d}   {j_E2sym[i]:7d}')

    C_AO_by_E = psi4.core.Matrix.array_interface(wfn.Ca_subset('AO','ALL'))[0]
    C_AO_by_sym = C_AO_by_E[:,j_E2sym]
    n_irrep = wfn.nirrep()
    n_mo_pi = wfn.nmopi()
    Ns = n_mo_pi.to_tuple()
    irrep_lst = []
    offset = 0
    for i in range(n_irrep):
        dim = Ns[i]
        irrep_lst.append(C_AO_by_sym[:,offset:offset+dim])
        offset += dim
    
    #C_AO_pi = psi4.core.Matrix.from_array(irrep_lst)
    return psi4.core.Matrix.from_array(irrep_lst)


def AO_to_SO(C_AO_pi, wfn):
    """
    transforms n_irrep blocks of Psi4-MOs from the Psi4-AO 
    to the Psi4-SO representation

    Parameters
    ----------
    C_AO_pi : psi4.core.Matrix C[AO:MO[irrep]]
    wfn : psi4.core.wavefunction (tested for RHF)

    Returns
    -------
    psi4.coreMatrix C[SO:MO[irrep]]

    """
    Ls = wfn.aotoso()
    n_irrep = wfn.nirrep()
    irrep_lst = []
    for sym in range(n_irrep):
        L, C = Ls.nph[sym], C_AO_pi.nph[sym]
        irrep_lst.append( np.matmul(np.transpose(L), C) )
    return psi4.core.Matrix.from_array(irrep_lst)


"""

This function doesn't change the outcome of the transformation.
We keep Psi4 order, it is simply a sum over AOs, no matter in which
order.

In other words, Cfour really uses a different SO order, and I have no clue
what that order would be.




"""


def AO_to_Cfour_SO(C_AO_pi, wfn, map_p2c, scale):
    """
    transforms n_irrep blocks of Psi4-MOs from the Psi4-AO to the 
    Cfour-SO representation

    Parameters
    ----------
    C_AO_pi : psi4.core.Matrix C[AO:MO[irrep]]
    wfn : psi4.core.wavefunction (tested for RHF)
    map_p2c: np.array[AO]
    scale:np.array[AO]
    Cfour4_MO[map_p2c[i]] = Psi4_MO[i]/scale[i]
    Cfour_AO_to_SO[map_p2c[i]] = Psi4_AO_to_SO[i]

    Returns
    -------
    psi4.coreMatrix C[SO:MO[irrep]]
    same dimensions as Psi4 MOs, but reordered and scaled

    """
    Ls = wfn.aotoso()
    n_irrep = wfn.nirrep()
    irrep_lst = []
    for sym in range(n_irrep):
        #L = Ls.nph[sym]
        L = psi4_to_c4(Ls.nph[sym], map_p2c, scale, use_scale=False)
        #C = C_AO_pi.nph[sym]
        C = psi4_to_c4(C_AO_pi.nph[sym], map_p2c, scale, use_scale=True)
        irrep_lst.append( np.matmul(np.transpose(L), C) )
    return psi4.core.Matrix.from_array(irrep_lst)


def psi4_to_c4(Cp4, map_p2c, scale, use_scale=True):
    """
    reorder the Psi4_Matrix C into Cfour order based on map_p2c[]
    Cfour4_C[map_p2c[i]] = Psi4_C[i]/scale[i]  [default, for MOs]
    Cfour4_C[map_p2c[i]] = Psi4_C[i]           [for AO_to_SO matrix]
    
    Parameters
    ----------
    Cp4 : np.array(nao, nmo) Psi4-matrix C    
    map_p2c : np.array(nao) mapping vector Psi4 to Cfour
    scale : np.array(nao) scaling vector Psi4 to Cfour
    use_scale: bool; don't use scale if reordering a AO-to-SO matrix

    Returns
    -------
    Cc4:  np.array(nao, nmo) Cfour4-matrix C

    """
    nbf = len(map_p2c)
    n_ao, n_mo = Cp4.shape
    if n_ao != nbf:
        msg = 'Error in psi4_to_c4: inconsistent AO-dimensions %d %d' % (nbf, n_ao)
        sys.exit(msg)
    
    Cc4=np.zeros((n_ao,n_mo))
    if use_scale:
        for i in range(n_ao):
            Cc4[map_p2c[i],:] = Cp4[i,:]/scale[i]
    else:
        for i in range(n_ao):
            Cc4[map_p2c[i],:] = Cp4[i,:]        
    return Cc4


def basis_mapping(basisset, verbose=1):
    """
    Computes the arrays needed for the operation:
    Cfour4_MO[map_p2c[i]] = Psi4_MO[i]/scale[i]    
    map_p2c = where to put MO coefficient i in the Cfour vector
    scale = how to scale MO coefficient i
    
    for Cfour, we assume all s, then all p, ...
    
    therefore in a Cfour_MO vector: s, s, s, px, px, py, py, pz, pz, d ...
    for every atom
    
    so for every atom we need first the offsets: shell_offset[l,m]
    then we can go over the atom-basis again and compute map_p2c[i]
    
    Parameters
    ----------
    basisset : psi4.core.basisset, Psi4 class for basis sets
    verbose : int, output level.

    Returns
    -------
    two vectors to be used in psi4_to_c4()
    map_p2c : int np.array; mapping Psi4 to Cfour (see above)
    scale : np.array; scaling factors Psi4 to Cfour (see above)
    """
    max_l = 5
    #max_m = 2*max_l + 1
    
    mol = basisset.molecule()
    n_atoms = mol.natom()
    nbf_atom = np.zeros(n_atoms, int)
    nbf = basisset.nbf()
    map_p2c = np.zeros(nbf, int)
    scale = np.zeros(nbf)

    """ 
      map from Psi4 to Cfour m-order  (trial and error)
      m_map[l,m], s_map[l,m]
      within a shell with fixed l:
      m_map = where to put MO coefficient m in the Cfour vector
      s_map = how to scale MO coefficient m
      p:    (0,1,2)         -> (2,0,1)   
      Psi4 uses (z,x,y) so 0(z) to 2, 1(x) to 0, 2(y) to 1 = (x,y,z)
      d:    (0,1,2,3,4)     -> (0,2,4,3,1)
      f:    (0,1,2,3,4,5,6) -> (2,0,1,6,4,3,5)
      and MO coefficient scaling (even more trial and error)
    """
    m_map = np.array(
        [[0, 0, 0, 0, 0, 0, 0],
         [2, 0, 1, 0, 0, 0, 0],
         [0, 2, 4, 3, 1, 0, 0],
         [2, 0, 1, 6, 4, 3, 5]]
        )
    sq12, sq24, sq40, sq60 = np.sqrt([12, 24, 40, 60])
    s_map = np.array(
        [[1, 0, 0, 0, 0, 0, 0],
         [1, 1, 1, 0, 0, 0, 0],
         [sq12, 1, 1, 2, 1, 0, 0],
         [sq60, sq40, sq40, 2, 1, sq24, sq24]]
        )                                      

    atom_offset = 0  # identical for Psi4_MO and Cfour_MO
    i_mo = 0

    for j_atom in range(mol.natom()):
        if verbose > 0:
            print(f'Atom {j_atom}')
        """ counter for shells for each angular momentum """
        nl = np.zeros(max_l, int)  
        nm = 2*np.arange(max_l)+1  # assume all shells are spherical

        """ initial sweep: establish shell_offset[l,m] """
        for i_shell in range(basisset.nshell_on_center(j_atom)):
            shell = basisset.shell(j_atom, i_shell)
            am = shell.amchar
            l = shell.am
            if verbose > 0:
                print(f"  Shell {i_shell}  {am} {l}")
            nl[l] += 1
        nbf_atom[j_atom] = np.sum(nl*nm)
        
        """ offset in the Cfour_MO vector: shell_offset[l,m] IMHO """
        shell_offset = np.zeros((6,9), int)
        rel_offset = 0
        shell_offset[0,0] = 0
        for l in range(1, max_l):
            rel_offset += nl[l-1]
            shell_offset[l,0] = rel_offset
            for m in range(1, 2*l+1):
                rel_offset += nl[l]
                shell_offset[l,m] = rel_offset
        if verbose > 1:
            print('shell_offset[l,m]:')
            print(shell_offset)

        """ second sweep: """
        nl = np.zeros(max_l, int)
        for i_shell in range(basisset.nshell_on_center(j_atom)):
            shell = basisset.shell(j_atom, i_shell)
            l = shell.am
            if verbose > 0:
                print(f'l = {l}, nl[l]={nl[l]}')
            for m in range(2*l+1):
                m_cfour = m_map[l,m]
                ptr = atom_offset + shell_offset[l,m_cfour] + nl[l]
                if verbose > 1:
                    print(f'  m={m} -> m={m_cfour}   ptr={ptr}')
                map_p2c[i_mo] = ptr
                scale[i_mo] = s_map[l,m]
                i_mo += 1
            nl[l] += 1
        
        atom_offset += nbf_atom[j_atom]

    if verbose > 0:    
        print(f'nbf = {nbf}   sum(nbf_atom)={np.sum(nbf_atom)}')
    
    return map_p2c, scale


def ao_offset(basisset, verbose=0):
    """
    nice, but works only for STO-3G
    
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
    

def Cfour_irrep_order(n_psi, group, verbose=1):
    """
    reorders the the vector n_per_irepp from Psi4 into Cfour order
    
    cfour[map[i]] = psi[i]
    
    Parameters
    ----------
    n_psi: np.array(int); number of orbitals per irrep in Psi4
    group : str; Schoenflies symbol; lower case
    verbose : output level

    Returns
    -------
    n_cfour; the reordered array according to Cfour convention
    
    --------
    C1, C2, Cs, Ci : nothing to do
    
    C2v: A1, A2, B1, B2 -> A1, B1, B2, A2
    C2h:   -> Ag, Au, Bu, Bg
    D2:   -> A, B2, B1, B3
    D2h:   -> Ag, B2u, B3u, B1g,   B1u, B2g, B3g, Au
    """
    #c2v = (0, 3, 1, 2)
    if len(n_psi) < 4:
        return n_psi


def read_oldmos(fname, nmos, RHF=True, verbose=1):
    """
    read OLDMOS to compare with created PSIMOS
    
    The number of MOs per irrep are known.
    We assume square matrices.

    Parameters
    ----------
    fname : str; file name
    nmos : np.array(int); number of MOs per irrep in Cfour-irrep order  
    RHF : Bool; False = UHF
    verbose : int; output level 

    Returns
    -------
    a list with array of MO coefficients in Cfour-irrep order
    """
    
    """ Compute number of expected lines """
    nirrep = len(nmos)
    lexp = 0
    for irrep in range(nirrep):
        lexp += nmos[irrep] * (nmos[irrep]//4 + 1)

    if verbose > 0:
        print('reading orbitals from '+fname)
    file = open(fname)
    lines = file.readlines()
    file.close()
    if len(lines) != lexp:
        print(f'{len(lines)} lines read from {fname}.')
        print(f'But {lexp} lines needed for n_mos_per_irrep vector.')
        sys.exit('Error')

    l = 0 # line counter
    Cs = []

    for h in range(nirrep):
        nc = nmos[h]
        #print(h,nc)
        Cs.append(np.zeros([nc,nc]))
        for jmo in range(0, nc, 4):
            #print(f'  jmo={jmo}')
            for ao in range(0, nc):
                #print(f'     jao={ao}, l={l}')
                words = lines[l].split()
                for i in range(len(words)):
                    Cs[h][ao,jmo+i] = float(words[i])
                l += 1

    return Cs


def read_oldmos_C1(fname, verbose = 1):
    """
    read a Cfour OLDMOS file
    ASSUMES nAOs = nMOs and tries to work out nAOs

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



def write_oldmos(fname, Cas, Cbs=None):
    """
    Write MOs in Cfour OLDMOS format
    That means packages of four MOs   
    C[0,0] C[0,1] C[0,2] C[0,3]
    C[1,0] C[1,1] C[1,2] C[1,3]
    C[2,0] C[2,1] C[2,2] C[2,3]
    ...    ...     ...    ...

    Format for each individual coefficient: 30.20E

    Write alpha and beta orbitals, or the alpha set twice.
    Writting the alpha set (RHF) twice doesn't affect a RHF guess
    and can serve as a UHF guess.

    Parameters
    ----------
    fname : str filename
    Cs : np.array with MO coefficients Cs[ao,mo]

    Returns
    -------
    None.
    """

    f = open(fname, 'w')
    nbf = Cas.shape[0]

    for j in range(0, nbf, 4):
        """ Normally we write groups of 4 MOs. The last group may be smaller. """
        ngr = 4
        if j + ngr >= nbf:
            ngr = nbf - j
        for ao in range(nbf):
            line = ''
            for igr in range(ngr):
                line += f"{Cas[ao,j+igr]:30.20E}"
            line += "\n"
            f.write(line)
    
    if Cbs is None:
        Cbs = Cas
    for j in range(0, nbf, 4):
        """ Normally we write groups of 4 MOs. The last group may be smaller. """
        ngr = 4
        if j + ngr >= nbf:
            ngr = nbf - j
        for ao in range(nbf):
            line = ''
            for igr in range(ngr):
                line += f"{Cbs[ao,j+igr]:30.20E}"
            line += "\n"
            f.write(line)
        
    f.close()
    return