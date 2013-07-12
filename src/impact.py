import trep
import numpy as np
import scipy
import scipy.optimize
from detect import *
from util import *


def linear_interpolate_impact(mvi):
    """ Finds a time and configuration when phi=0 based on linear interpolation.
        this is used to get a guess to seed the root-solver. """

    set_system_state(mvi, 1)
    phi0 = phi(mvi)
    set_system_state(mvi, 2)
    phi1 = phi(mvi)

    ti = mvi.t1 - (phi0/(phi1-phi0))*(mvi.t2-mvi.t1)
    frac = (ti-mvi.t1)/(mvi.t2-mvi.t1)
    qi = frac*(mvi.q2-mvi.q1)+mvi.q1

    return ti, qi


def solve_impact_time(mvi):
    """ Finds the time and configuration at impact. """

    n, m, nk = mvi.nd, mvi.nc, mvi.nk

    # Linearly interpolate the kinematic configuration at the moment of impact.
    qk1, qk2 = mvi.q1[n:], mvi.q2[n:]
    def kin_config(t):
        return (t-mvi.t1)/(mvi.t_end-mvi.t1)*(qk2-qk1)+qk1

    set_system_state(mvi, 1)
    Dh1  = Dh(mvi)
    Dh1T = np.transpose(Dh1)

    def func(x):
                
        mvi.t2 = x[-1]
        mvi.q2 = np.append(x[:n], kin_config(x[-1]))
        mvi.lambda1 = x[(n+nk):(n+nk+m)]

        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T[:n],mvi.lambda1)

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f3 = phi(mvi)
        return np.hstack((f1, f2, f3))

    t_guess, q_guess = linear_interpolate_impact(mvi)
    x0 = np.hstack((q_guess[:n], mvi.lambda1, t_guess))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True)

    if ier != 1:
        # The root solver will fail for exactly simultaneous collisions since setting
        # t1 = t2 will return nan for the system's velocity. Can't handle this right now, 
        # but can at least check for it.
        if mvi.t1 == x[-1]:
            print 'Simultaneous impacts detected.'
        raise trep.ConvergenceError("Find impact failed to converge.")      

    mvi.q2 = np.append(x[:n], kin_config(x[-1]))
    mvi.t2 = x[-1]
    mvi.lambda1 = x[(n+nk):(n+nk+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lam1': mvi.lambda1, 'lam0': mvi.lambda0}


def plastic_impact(mvi):
    """ Finds the time and configuration at impact. """

    n, m, nk = mvi.nd, mvi.nc, mvi.nk

    mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    mvi.set_single_state(1, mvi.t2, mvi.q2, mvi.p2)
    mvi.t2 = mvi.t_end
    mvi.lambda0 = mvi.lambda1

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))
    Dphi1 = Dphi(mvi)
    D2h1  = D2h(mvi) 

    mvi.set_midpoint()
    tau = D4L2(mvi)

    def func(x):
        mvi.q2 = np.append(x[:n], mvi.qk_end)
        mvi.lambda1 = x[(n+nk):(n+nk+m)]
        lambda_c = x[-2]
        Ec = x[-1]

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f4 = phi(mvi)

        mvi.set_midpoint()
        f1 = (mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T[:n],mvi.lambda1)
              -np.transpose(Dphi1)*lambda_c)[0]
        f3 =  (tau+D3L2(mvi)-np.dot(D2h1,mvi.lambda1)+Ec)
        return np.hstack((f1, f2, f3, f4))

    x0 = np.hstack((mvi.q2[:n], mvi.lambda1, 0.0, 0.0))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True)

    if ier != 1:
        raise trep.ConvergenceError("Impact update failed to converge.")

    mvi.q2 = np.append(x[:n],mvi.qk_end)
    mvi.lambda1 = x[(n+nk):(n+nk+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1, 'lambdac': x[-2], 'Ec': x[-1]}
    

def solve_impact(mvi, impact_frames):
    """ When given one or more impacts, finds the first impact to occur and the state
        at that point. Returns the frame and state. """

    q2, p2, t2, lam1 = mvi.q2, mvi.p2, mvi.t2, mvi.lambda1

    impacts = {}
    for frame in impact_frames:
        mvi.surface.set_frame(frame)
        s = solve_impact_time(mvi)
        impacts.update({frame.name: s})
        mvi.q2, mvi.p2, mvi.t2, mvi.lambda1 = q2, p2, t2, lam1

    frame_name = impacts.keys()[np.argmin([impacts[k]['t'] for k in impacts.keys()])]
    frame = mvi.system.get_frame(frame_name)
    impact = impacts[frame_name]
    
    return frame, impact    


def impact_update(mvi, impact_frame, t, q, p, lam0, lam1):
    """ Given a state at the time of impact, this function solves the impact map and 
        adds the new constraint to the system. """  
     
    # Set system config to impact config. Solve impact update.
    mvi.t2, mvi.q2, mvi.p2, mvi.lambda0, mvi.lambda1 = t, q, p, lam0, lam1
    mvi.surface.set_frame(impact_frame)
    next_state = plastic_impact(mvi)
    state = save_state(mvi)

    # Add the constraint to the system.
    mvi.constrained_frames.append(impact_frame)
    mvi.surface.add_constraint_to_system(impact_frame)

    # Adding constraints can sometimes clear the inteagrator's states. Here, we restore the
    #   current state. 
    restore_state(mvi, state)
    mvi.lambda0 = np.append(mvi.lambda0, [0.0]) 
    mvi.lambda1 = np.append(next_state['lambda1'], next_state['lambdac'])           
