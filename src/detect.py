from util import *

def detect_impacts(mvi):
    """ Detects which surfaces, if any, have been crossed. """
    set_system_state(mvi, 2)
    impact_list = []
    for surface in mvi.inactive_surfaces: 
        if surface.phi() < 0:
            impact_list.append(surface)     
    return impact_list

def detect_releases(mvi):
    if mvi._releases_off:
        return []
    release_list = []
    lam2 = mvi.lambda2c
    lam1 = mvi.lambda1c
    for surface in mvi.active_surfaces:
        if lam2[surface.index]*surface.sign < 0.0:
            release_list.append(surface)
    return release_list  

def sign_change(a, b):
    if (a > 0) and (b < 0):
        return True
    if (a < 0) and (b > 0):
        return True
    return False

'''
def detect_releases2(mvi):
    """ Detects which nodes, if any, should be released. This is based
        on a sign change in the force of constraint. """ 
    if mvi._releases_off:
        return []
    release_list = []
    for surface in mvi.active_surfaces:
        con_index = mvi.system.get_constraint(surface).index
        if sign_change(mvi, con_index):
            release_list.append(surface)           
    return release_list

def sign_change(mvi, i):
    """ Returns True if the i-th constraint force has changed sign, False otherwise. """
    if (mvi.lambda0[i] > 0) and (mvi.lambda1[i] < 0):
        return True
    if (mvi.lambda0[i] < 0) and (mvi.lambda1[i] > 0):
        return True
    return False     
'''    
