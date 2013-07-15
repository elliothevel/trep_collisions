import trep
import numpy as np
import trep_collisions as tc


class PuppetString:
    def __init__(self, constraint):
        self.system = constraint.system
        self.constraint = constraint
        self.active = True
      
    def activate(self):
        self.system._constraints += (self.constraint,)
        self.active = True 
   
    def deactivate(self):
        new_constraints = tuple()
        for con in self.system.constraints:
            if con != self.constraint: 
                new_constraints += (con, )
        self.system._constraints = new_constraints
        self.active = False
           
                                
def detect_string_releases(pmvi):
    """ Checks if a string constraint which is currently active should
        be released because of a change in sign of the constraint multiplier. """    
    release_strings = []
    for string in pmvi.puppet_strings:
        if string.active:
            if tc.detect.sign_change(mvi, string.constraint.index):
                release_strings.append(string)
    return release_strings


def check_string_tension(pmvi):
    """ Checks if a string constraint which is currently not active should
        be activated because the distance has exceeded the distance specified 
        by the constraint. """    
    catch_strings = []
    for string in pmvi.puppet_strings:
        if not string.active:
            c = string.constraint
            if c.get_actual_distance() > c.distance:
                catch_strings.append(string)
    return catch_strings            


def get_catch_time(mvi, catch_strings):
    """ Determines the time at which one or more string constraints should be
        activated as calculated by linear interpolation. """

    d1, d2 = [], []

    # Set the system to the first state and find the actual distances.
    tc.util.set_system_state(mvi, 1)
    for string in catch_strings:
        d1.append(string.constraint.get_actual_distance()-string.constraint.distance)

    # Set the system to the second state and find the actual distances.
    tc.util.set_system_state(mvi, 2)
    for string in catch_strings:   
        d2.append(string.constraint.get_actual_distance()-string.constraint.distance)
        
    # Linearly interpolate to find the catching time.
    times = [mvi.t1 - d1[i]/(d2[i]-d1[i])*(mvi.t2-mvi.t1) for i in range(len(d1))]
    return np.mean(times) 
    
    
def solve_string_catch(mvi, catch_strings, t):

    # Reset the integrator to state 1.
    qk1, qk2 = mvi.q1[mvi.nd:], mvi.q2[mvi.nd:]
    qk = (t-mvi.t1)/(mvi.t2-mvi.t1)*(qk2-qk1) + qk1
    mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
    mvi.set_single_state(1, mvi.t0, mvi.q0, mvi.p0)
    mvi.set_single_state(0, mvi.t0, mvi.q0, mvi.p0)
    mvi.lambda1 = mvi.lambda0     

    # Step to the time of the string catch.
    mvi.step(t, k2=tuple(qk))
    state = tc.util.save_state(mvi)

    # Activate the strings.
    mvi.system.hold_structure_changes()
    activate_strings(catch_strings)
    mvi.system.resume_structure_changes()

    # Restore state and step to the final time.
    tc.util.restore_state(mvi, state)
    mvi.lambda0 = np.append(state['lambda1'], np.zeros(len(catch_strings)))
    mvi.lambda1 = np.append(state['lambda1'], np.zeros(len(catch_strings)))
    mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    mvi.step(mvi.t_end, k2=tuple(qk))
     
                                       
def deactivate_strings(string_list):
    for string in string_list:
        string.deactivate()


def activate_strings(string_list):
    for string in string_list:
        string.activate()    


def release_update_with_strings(mvi, frames, strings, t, q, p, lam0):
    """ Removes regular constraints and string constraints. """

    indices = ([mvi.system.get_constraint('constraint_' + frame.name).index for frame in frames]
               +[string.constraint.index for string in strings])

    # Save the last two states of the integrator.
    t2, q2, p2 = t, q, p
    t1, q1, p1 = mvi.t0, mvi.q0, mvi.p0
    lam1 = np.delete(lam0, indices)

    # Remove the constraints from the system and go back to mvi.t1.
    mvi.system.hold_structure_changes()
    tc.release.remove_constraints(mvi, frames)
    deactivate_strings(strings)
    mvi.system.resume_structure_changes()
    
    # Add the state back.
    mvi.set_single_state(2, t2, q2, p2)
    mvi.set_single_state(1, t1, q1, p1)
    mvi.set_single_state(0, t1, q1, p1)
    mvi.lambda0, mvi.lambda1 = lam1, lam1

    # Remove the constraint and simulate to the end of the current timestep. 
    mvi.constrained_frames = [f for f in mvi.constrained_frames if (f not in frames)]
    if mvi.t_end != mvi.t2:
        mvi.step(mvi.t_end, k2=q2[mvi.nd:])             
