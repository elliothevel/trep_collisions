import trep
import numpy as np
import scipy
import scipy.optimize
from detect import *
from util import *
from impact import *


LAMBDA_TOL = 1e-5


def solve_release(mvi, release_surfaces, method):
    """ Solves for the release times of several frames. """

    # Endpoint method: release occurs at state 1.
    if method == 'endpoint':
        state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lambda1': mvi.lambda0}
        return release_surfaces, state
    
    # Interp method: linearly interpolate the time of release.
    elif method == 'interp':

        times = []
        for surf in release_surfaces:
            times.append(linear_interpolate_release(mvi, surf.index)[0])
        tr = np.mean(times)

        if abs(mvi.t1-tr)<1e-5:
            mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
            mvi.lambda1 = mvi.lambda0
            state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}
            return release_surfaces, state

        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.step(tr, k2=mvi.kin_config(tr))
        state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}
        return release_surfaces, state
    
    # Root method: root-solve for time of release. This finds the first time of
    # release among all release surfaces, then checks for multiple releases.
    elif method == 'root':
        q2, p2, t2, lam1 = mvi.q2, mvi.p2, mvi.t2, mvi.lambda1
        releases = {}
        for surface in release_surfaces:
            s = root_solve_release(mvi, surface)  
            releases.update({surface: s})
            mvi.q2, mvi.p2, mvi.t2, mvi.lambda1 = q2, p2, t2, lam1 
        
        release_surface = releases.keys()[np.argmin([releases[k]['t'] for k in releases.keys()])]
        release_state = releases[surface]

        # Check for multiple releases.
        surfaces = [release_surface]
        mvi.set_single_state(2,  release_state['t'],  release_state['q'],  release_state['p'])
        mvi.lambda1 = release_state['lambda1']
        set_system_state(mvi, 2)
        for surface in release_surfaces:
            if surface != release_surface:
                if abs(mvi.lambda2c[surface.index]) < LAMBDA_TOL:
                    surfaces.append(surface)  
        
        return surfaces, release_state 
    

def linear_interpolate_release(mvi, j):
    """ Linearly interpolates from constraint forces the time and configuration at release. """

    lam1 = mvi.lambda1c[j]
    lam2 = mvi.lambda2c[j]

    tr = mvi.t1 - (lam1/(lam2-lam1))*(mvi.t2-mvi.t1)
    frac = (tr-mvi.t1)/(mvi.t2-mvi.t1)
    qr = frac*(mvi.q2-mvi.q1)+mvi.q1

    return tr, qr
 

def root_solve_release(mvi, con):
    """ Solves for the state of the system at the time of release of constraint con. 
        This is when the continuous-time Lagrange multiplier is zero.

        This modifies lambda1 and state 2.
    """
    nd, m = mvi.nd, mvi.nc

    # If the first state is already close enough to release, step integrator back
    # and use the state as the release state.
    if abs(mvi.lambda1c[con.index]) < LAMBDA_TOL:
        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.lambda1 = mvi.lambda0
        return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))  

    # The root-solving problem is to find (DEL, h(qr, tr), lambda_{cont}) = 0. 
    def func(x):
        mvi.t2 = x[-1]
        mvi.lambda1 = x[nd:-1]
        mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))

        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f3 = mvi.system.lambda_(constraint=con)

        return np.hstack((f1, f2, f3))


    def fprime(x):
        mvi.t2 = x[-1]
        mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))
        mvi.lambda1 = x[(nd):(nd+m)]
        
        mvi.set_midpoint()
        Df11 = D2D1L2(mvi) + D2fm2(mvi)
        Df13 = D4D1L2(mvi) + D4fm2(mvi)
        
        set_system_state(mvi, 2)
        Df21 = Dh(mvi)
        Df23 = D2h(mvi).reshape((mvi.nc, 1))
        lam_dq = (mvi.system.lambda_dq(constraint=con)+
                 mvi.system.lambda_ddq(constraint=con)/(mvi.t2-mvi.t1))[:nd]
        lam_dt = -1.0/(mvi.t2-mvi.t1)*np.dot(mvi.system.dq[:nd], 
                       mvi.system.lambda_ddq(constraint=con)[:nd])

        Df1 = np.hstack((Df11, -Dh1T, Df13))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc)), Df23))
        Df3 = np.hstack((lam_dq, np.zeros(m), lam_dt))

        return np.vstack((Df1, Df2, Df3))
   

    t_guess, q_guess = linear_interpolate_release(mvi, con.index)
    x0 = np.hstack((q_guess[:nd], mvi.lambda0, t_guess))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Release solver failed to converge.") 
  
    mvi.set_midpoint()
    mvi.p2 = D2L2(mvi)
    mvi.q2 = x[:nd]
    mvi.t2 = x[-1]
    mvi.lambda1 = x[nd:-1]

    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}


def plastic_impact_to_release(mvi, impact_surfaces, con):
    pass


def remove_constraints(mvi, surfaces):
    """ Removes the constraints on frames from the system. """
    mvi.system.hold_structure_changes()
    for surface in surfaces:
        surface.deactivate_constraint()
    mvi.system.resume_structure_changes()    


def release_update(mvi, surfaces, release_state, stepafter=True):
    """ Given a state at the point of release, this function removes the constraint(s) from the
        system and intagrates to the next time step. """ 
    
    # Save the last two states of the integrator.
    t2, q2, p2 = release_state['t'], release_state['q'], release_state['p']
    t1, q1, p1 = mvi.t1, mvi.q1, mvi.p1
    indices = [surface.index for surface in surfaces]
    lam0 = np.delete(mvi.lambda0, indices)
    lam1 = np.delete(release_state['lambda1'], indices)

    # Remove the constraints from the system and go back to mvi.t1.
    remove_constraints(mvi, surfaces)
    
    # Add the state back.
    mvi.set_single_state(2, t2, q2, p2)
    mvi.set_single_state(1, t1, q1, p1)
    #mvi.set_single_state(0, t1, q1, p1)
    mvi.lambda0, mvi.lambda1 = lam0, lam1

    # Remove the constraint and simulate to the end of the current timestep. 
    if stepafter:
        if mvi.t_end != mvi.t2:
            mvi.lambda0 = mvi.lambda1
            mvi.step(mvi.t_end, k2=mvi.kin_config(mvi.t_end))    


def contact_transfer_update(mvi, impact_surfaces, release_surfaces, transfer_state):
    """ This solves the problem of constraints being added and removed at the same time. """

        # Save the last two states of the integrator.
    t2, q2, p2 = release_state['t'], release_state['q'], release_state['p']
    t1, q1, p1 = mvi.t1, mvi.q1, mvi.p1
    indices = [surface.index for surface in surfaces]
    lam0 = np.delete(mvi.lambda0, indices)
    lam1 = np.delete(release_state['lambda1'], indices)

    # Remove the constraints from the system and go back to mvi.t1.
    remove_constraints(mvi, surfaces)
    
    # Add the state back.
    mvi.set_single_state(2, t2, q2, p2)
    mvi.set_single_state(1, t1, q1, p1)
    #mvi.set_single_state(0, t1, q1, p1)
    mvi.lambda0, mvi.lambda1 = lam0, lam1

    # Remove the constraint and simulate to the end of the current timestep. 
    if mvi.t_end != mvi.t2:
        mvi.lambda0 = mvi.lambda1
        mvi.step(mvi.t_end, k2=mvi.kin_config(mvi.t_end)) 

    impact_update(mvi, impact_surfaces, transfer_state)    

        
    
