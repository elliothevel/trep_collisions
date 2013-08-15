import trep
import numpy as np
import scipy
import scipy.optimize
from detect import *
from util import *
from impact import *


LAMBDA_TOL = 1e-5
TIME_TOL   = 1e-5


def find_first_release(mvi, release_surfaces, method):
    """ Solves for the release times of several frames. This method does not
        modify the states of the integrator.
    """

    # Endpoint method: release occurs at state 1.
    if method == 'endpoint':
        state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lambda1': mvi.lambda0}
        return release_surfaces, state
    
    # Interp method: linearly interpolate the time of release.
    elif method == 'interp':

        times = []
        for surface in release_surfaces:
            times.append(linear_interpolate_release(mvi, surface.index)[0])
        tr = np.mean(times)

        ## if abs(mvi.t1-tr)<TIME_TOL:
        ##    state = {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lambda1': mvi.lambda0}
        ##    return release_surfaces, state

        initial_state = save_state(mvi)
        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)

        if not mvi._state1_impacts:
            mvi.step(tr, k2=mvi.kin_config(tr))
        else:
            tend = mvi.t_end
            mvi.t_end = tr

            if abs(tr-mvi.t1)<TIME_TOL:
                transfer_state = {'p': mvi.p1, 't':mvi.t1, 'q': mvi.q1, 'tau': mvi.tau2, 'lambda1': mvi.lambda1}
                constraint_transition_update(mvi, release_surfaces, mvi._state1_impacts, transfer_state)
            else:    
                plastic_impact(mvi)
            mvi.t_end = tend

        state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}
        restore_state(mvi, initial_state)

        return release_surfaces, state
    
    # Root method: root-solve for time of release. This finds the first time of
    # release among all release surfaces, then checks for multiple releases.
    elif method == 'root':
        state = save_state(mvi)
        releases = {}
        for surface in release_surfaces:
            s = solve_release_time(mvi, surface)  
            releases.update({surface: s})
            restore_state(mvi, state)

        release_surface = releases.keys()[np.argmin([releases[k]['t'] for k in releases.keys()])]
        release_state = releases[surface]

        # Check for multiple releases.
        surfaces = [release_surface]
        mvi.set_single_state(2,  release_state['t'],  release_state['q'],  release_state['p'])
        mvi.lambda1 = release_state['lambda1']
        set_system_state(mvi, 2)
        mvi.system.set_state()
        for surface in release_surfaces:
            if surface != release_surface:
                if abs(mvi.lambda2c[surface.index]) < LAMBDA_TOL:
                    surfaces.append(surface)  
        restore_state(mvi, state)            
        
        return surfaces, release_state 
    

def linear_interpolate_release(mvi, j):
    """ Linearly interpolates from constraint forces the time and configuration at release. """

    lam1 = mvi.lambda1c[j]
    lam2 = mvi.lambda2c[j]

    if (lam1 < 0) and (lam2 < 0):
        print lam1, lam2 
        raise Exception("Bad release interpolation.")

    tr = mvi.t1 - (lam1/(lam2-lam1))*(mvi.t2-mvi.t1)
    frac = (tr-mvi.t1)/(mvi.t2-mvi.t1)
    qr = frac*(mvi.q2-mvi.q1)+mvi.q1

    return tr, qr
 

def solve_release_time(mvi, con):
    """ Solves for the state of the system at the time of release of constraint con. 
        This is when the continuous-time Lagrange multiplier is zero.

        This method modifies lambda1 and state 2.
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
    D2h1 = D2h(mvi)

    t_guess, q_guess = linear_interpolate_release(mvi, con.index)

    # No impacts at state 1 - solve normal release time problem.
    if not mvi._state1_impacts:

        x0 = np.hstack((q_guess[:nd], mvi.lambda0, t_guess))

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

    # Impacts have occurred at state 1 - solve impact update with additional unknown 
    # release time.
    else:

        print '\ttesting'
        impact_surfaces = mvi._state1_impacts
        p = len(impact_surfaces)

        tau = mvi.tau2      # Is set by most recent impact solver.

        x0 = np.hstack((q_guess[:nd], mvi.lambda1, t_guess, np.dot(D2h1, mvi.lambda1)-tau))

        def func(x):  
            mvi.t2 = x[nd+m]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            Ec = x[-1]

            mvi.set_midpoint()
            f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)
            f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

            set_system_state(mvi, 2)
            f2 = h(mvi)
            f4 = mvi.system.lambda_(constraint=con)
               
            return np.hstack((f1, f2, f3, f4))

        def fprime(x):
            mvi.t2 = x[nd+m]
            mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
            mvi.lambda1 = x[nd:(nd+m)]
            Ec = x[-1]
        
            mvi.set_midpoint()
            Df11 = D2D1L2(mvi) + D2fm2(mvi)
            Df13 = D4D1L2(mvi) + D4fm2(mvi)
            Df31 = D2D3L2(mvi)
            Df33 = D4D3L2(mvi)
        
            set_system_state(mvi, 2)
            Df21 = Dh(mvi)
            Df23 = D2h(mvi).reshape((m, 1))
            Df32 = -D2h(mvi)
            Df41 = D2LAM2(mvi, con)
            Df43 = D4LAM2(mvi, con)

            Df1 = np.hstack((Df11, -Dh1T, Df13, np.zeros((nd,1))))
            Df2 = np.hstack((Df21, np.zeros((m, m)), Df23, np.zeros((m,1))))
            Df3 = np.hstack((Df31, Df32, Df33, (1.0,)))
            Df4 = np.hstack((Df41, np.zeros(m), Df43, np.zeros(1)))

            return np.vstack((Df1, Df2, Df3, Df4))

   
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True)#, fprime=fprime)

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Release solver failed to converge.") 

    print 'solved', func(x)
    print x0[nd+m], x[nd+m]
    print mvi.t1, mvi.t_end
    mvi.t2 = x[nd+m]
    mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}


def remove_constraints(mvi, surfaces):
    """ Removes the constraints on frames from the system. """

    state = save_state(mvi)
    indices = [surface.index for surface in surfaces]
    lam0 = np.delete(mvi.lambda0, indices)
    lam1 = np.delete(mvi.lambda1, indices)

    mvi.system.hold_structure_changes()
    for surface in surfaces:
        surface.deactivate_constraint()
    mvi.system.resume_structure_changes() 

    restore_state(mvi, state)
    mvi.lambda0 = lam0
    mvi.lambda1 = lam1


def release_update(mvi, surfaces, release_state):
    """ Given a state at the point of release, this function removes the constraint(s) from the
        system and intagrates to the next time step. """ 
    
    # Save the last two states of the integrator.
    mvi.set_single_state(2, release_state['t'], release_state['q'], release_state['p'])
    t2, q2, p2 = release_state['t'], release_state['q'], release_state['p']

    # Remove the constraints from the system and go back to mvi.t1.
    remove_constraints(mvi, surfaces)
    
    # Remove the constraint and simulate to the end of the current timestep. 
    if mvi.t_end != mvi.t2:
        mvi.lambda0 = mvi.lambda1
        mvi._dq1 = Dq(mvi)
        mvi.step(mvi.t_end, k2=mvi.kin_config(mvi.t_end)) 
        mvi._state1_impacts = []


def constraint_transition_update(mvi, impact_surfaces, release_surfaces, transfer_state):
    """ This solves the problem of constraints being added and removed at the same time. """
    mvi.set_single_state(2, transfer_state['t'], transfer_state['q'], transfer_state['p']) 
    mvi.lambda1 = transfer_state['lambda1']
    mvi.set_midpoint()
    transfer_state['tau'] = D4L2(mvi)

    # Remove the constraints from the system.
    remove_constraints(mvi, release_surfaces)

    # Apply the impact map to the end of the interval.
    impact_update(mvi, impact_surfaces, transfer_state)    

        
    
