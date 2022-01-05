#!/usr/bin/env python3


"""
@author T. Sommerfeld 

Jan 5 2022

"""

def cfour_comp_sym_and_geo(fout):
    """
    extract from a Cfour output:
        the computational point group 
        the coordinates used in the calculation
        
    fout: filename of a Cfour output file
    
    returns: 
        str sym  (the point group)
        str coor (the computational geometry in Bohr)
    """
    
    lines = readfile(fout)
    
    l, n = 0, len(lines)
    
    coor = ''
    
    while l < n:
        if 'The computational point group is' in lines[l]:
            words = lines[l].split()
            sym = words[5]
            break
        l += 1
    
    while l < n:
        if 'Coordinates used in calculation (QCOMP)' in lines[l]:
            break
        l += 1
        
    l += 5
    while l < n:
        words = lines[l].split()
        if len(words) == 5:
            atom = (words[0], words[2], words[3], words[4])
            coor += '%-2s   %11s  %11s  %11s\n' % atom 
            l += 1
        else:
            break

    return sym, coor
    

def readfile(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    return lines







