# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 13:38:30 2022

@author: mlols
"""


def hk_modes(hier_num):
    """
    Generate modes in the HK hierarchy.

    Parameters
    ----------
    hier_num : int
        Number in the HK hierarchy (hier_num = n means the nth model).

    Returns
    -------
    p_modes : list
        List of psi modes, represented as tuples.
        Each tuple contains the horizontal and vertical wavenumbers.
    t_modes : list
        List of theta modes, represented as tuples.
        
    Examples
    --------
    >>> p_modes, t_modes = hk_modes(1)
    >>> print(p_modes)
    [(0, 1), (1, 1)]
    >>> print(t_modes)
    [(0, 2), (1,1)]
    >>> p_modes, t_modes = hk_modes(2)
    >>> print(p_modes)
    [(0, 1), (0, 3), (1, 1), (1, 2)]
    >>> print(t_modes)
    [(0, 2), (0, 4), (1,1), (1, 2)]
    """
    p_modes = [(0,1), (1,1)] #Base model
    t_modes = [(0,2), (1,1)]
    
    current_pair = (1,1) #
    
    #Add pairs of modes. Fills shells of constant L1 norm by increasing norm.
    #Within shell, fills modes in descending lexicographical order.
    #When starting new shell, add modes of zero horiz wavenumber.
    for i in range(1, hier_num):

        if current_pair[1] == 1:
            level = current_pair[0]+1
            current_pair = (1, level)
            p_modes.append((0, level*2-1))
            t_modes.append((0, level*2))
        else:
            current_pair = (current_pair[0]+1, current_pair[1]-1)
            
        p_modes.append(current_pair)
        t_modes.append(current_pair)
        
    p_modes.sort()
    t_modes.sort()
                
    return p_modes, t_modes