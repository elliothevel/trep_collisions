from util import *


def detect_impacts(mvi):
    """ Detects which nodes, if any, have made contact. """
    set_system_state(mvi, 2)
    frame_list = []
    for frame in mvi.impact_frames:
        if frame not in mvi.constrained_frames:
            mvi.surface.set_frame(frame)                                   
            if phi(mvi) < 0:
                frame_list.append(frame)                                   
    return frame_list
        
 
def detect_releases(mvi):
    """ Detects which nodes, if any, should be released. This is based
        on a sign change in the force of constraint. """    
    release_list   = []
    for frame in mvi.constrained_frames:
        con_index = mvi.system.get_constraint('constraint_'+frame.name).index
        if sign_change(mvi, con_index):
            release_list += [frame]     
    return release_list


def sign_change(mvi, i):
    """ Returns True if the i-th constraint force has changed sign, False otherwise. """
    if (mvi.lambda0[i] > 0) and (mvi.lambda1[i] < 0):
        return True
    if (mvi.lambda0[i] < 0) and (mvi.lambda1[i] > 0):
        return True
    return False       
