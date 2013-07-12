import trep
import numpy as np
import scipy
import scipy.optimize
from util import *


def solve_release(mvi, release_frames, root_solve):
    """ Solves for the release times of several frames. """

    if not root_solve:
        state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lam0': mvi.lambda0}
        return release_frames, state
    
    q1, p1, t1, lam0 = mvi.q1, mvi.p1, mvi.t1, mvi.lambda0
    releases = {}
    for frame in release_frames:
        j = mvi.system.get_constraint('constraint_' + frame.name).index 
        s = root_solve_release(mvi, j)  
        releases.update({frame.name: s})
        mvi.q1, mvi.p1, mvi.t1, mvi.lambda0 = q1, p1, t1, lam0 
        
    frame_name = releases.keys()[np.argmin([releases[k]['t'] for k in releases.keys()])]
    frame = mvi.system.get_frame(frame_name)
    release = releases[frame_name]
        
    return [frame], release 
    

def linear_interpolate_release(mvi, j):
    """ Linearly interpolates from constraint forces the time and configuration at release. """

    tr = mvi.t0 - (mvi.lambda0[j]/(mvi.lambda1[j]-mvi.lambda0[j]))*(mvi.t2-mvi.t0)
    frac = (tr-mvi.t0)/(mvi.t2-mvi.t0)
    qr = frac*(mvi.q2-mvi.q0)+mvi.q0

    return tr, qr
 

def root_solve_release(mvi, j):
    """ Solves for the state of the system at the time of release of constraint j. """

    n, m = mvi.nd, mvi.nc
    con = mvi.system.constraints[j]

    t1, p1, q1 = mvi.t1, mvi.p1, mvi.q1
    mvi.t1, mvi.p1, mvi.q1 = mvi.t0, mvi.p0, mvi.q0

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))

    # Get the right sign for the optimization. Since the constraints can change sign in
    #   either direction, the optimization will be either a minimization or a maximization. 
    if mvi.lambda1[j] >= 0:
        opt_sign = 1.0
    else:
        opt_sign = -1.0  

    # A built-in scipy solver is used. The goal is to maximize lambda_j while satisfying the dynamics.
    def func(x):
        return opt_sign*x[n+j]

    def fprime(x):
        fp = np.zeros(len(x))
        fp[n+j] = opt_sign*1.0             
        return fp

    def f_eqcons(x):
        mvi.q2, mvi.t2, mvi.lambda0 = x[:n], x[-1], x[n:(n+m)] 
        f = np.zeros(m+n)
        mvi.system.set_state(q=x[:n], t=x[-1])
        f2 = h(mvi)
        mvi.set_midpoint()
        f1 = (mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T[:n],mvi.lambda0))
        return np.hstack((f1, f2))

    bounds = ([(-10*abs(q),10*abs(q)) for q in mvi.q2[:n]] +
              [(-10*abs(l),10*abs(l)) for l in mvi.lambda1] + [(mvi.t0,mvi.t2)])    

    x0 = np.hstack((q1, mvi.lambda0, t1))
    (x_opt, fx, its, imode, smode) = scipy.optimize.fmin_slsqp(func, x0, f_eqcons=f_eqcons,
                                                               fprime=fprime, acc=mvi.tolerance, iter=1000,
                                                               bounds=bounds, iprint=0, full_output=True)

    if imode != 0:
        raise trep.ConvergenceError("Release solver failed to converge.")           
  
    mvi.set_midpoint()
    mvi.p1 = D2L2(mvi)
    mvi.q1 = x_opt[:n]
    mvi.t1 = x_opt[-1]
    mvi.lambda0 = x_opt[n:(n+m)]

    return {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lam0': mvi.lambda0}


def remove_constraints(mvi, frames):
    """ Removes the constraints on frames from the system. """
    constraints = [mvi.system.get_constraint('constraint_'+frame.name) for frame in frames]	

    new_constraints = tuple()
    for old_constraint in mvi.system.constraints:
        if old_constraint not in constraints:
            new_constraints += (old_constraint,)

    # Reset the list of constraints. Let trep know that the system
    #   has changed its structure (not sure if needed, but might help)        
    mvi.system._constraints = new_constraints
    mvi.system._structure_changed()


def remove_constraint_from_list(mvi, old_frame):
    """ Removes a released frame from the integrator's list of constrained frames. """
    mvi.constrained_frames = [f for f in mvi.constrained_frames if not (f == old_frame)]


def release_update(mvi, frames, t, q, p, lam0):
    """ Given a state at the point of release, this function removes the constraint from the
        system and intagrates to the next time step. """

    indices = [mvi.system.get_constraint('constraint_' + frame.name).index for frame in frames]

    # Save the last two states of the integrator.
    t2, q2, p2 = t, q, p
    t1, q1, p1 = mvi.t0, mvi.q0, mvi.p0
    lam1 = np.delete(lam0, indices)

    # Remove the constraints from the system and go back to mvi.t1.
    mvi.system.hold_structure_changes()
    remove_constraints(mvi, frames)
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
