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
    """ Finds the time and configuration at impact.
          when passed to this function, mvi should have following states:
            0 - two states before impact
            1 - state just before impact
            2 - invalid configuration to be corrected
    """

    nd, m = mvi.nd, mvi.nc

    # Linearly interpolate the kinematic configuration at the moment of impact.
    qk1, qk2 = mvi.q1[nd:], mvi.q2[nd:]
    def kin_config(t):
        return (t-mvi.t1)/(mvi.t_end-mvi.t1)*(qk2-qk1)+qk1

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))

    def func(x):
                
        mvi.t2 = x[-1]
        mvi.q2 = np.append(x[:nd], kin_config(x[-1]))
        mvi.lambda1 = x[nd:(nd+m)]

        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f3 = phi(mvi)

        return np.hstack((f1, f2, f3))

    def fprime(x):
        mvi.t2 = x[-1]
        mvi.q2 = np.append(x[:nd], kin_config(x[-1]))
        mvi.lambda1 = x[(nd):(nd+m)]
        
        mvi.set_midpoint()
        Df11 = D2D1L2(mvi) + D2fm2(mvi)
        Df13 = D4D1L2(mvi) + D4fm2(mvi)
        
        set_system_state(mvi, 2)
        Df21 = Dh(mvi)
        Df31 = Dphi(mvi)
        Df23 = D2h(mvi).reshape((mvi.nc, 1))

        Df1 = np.hstack((Df11, -Dh1T, Df13))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc)), Df23))
        Df3 = np.hstack((Df31, np.zeros(mvi.nc+1)))

        return np.vstack((Df1, Df2, Df3))
        
    t_guess, q_guess = linear_interpolate_impact(mvi)
    x0 = np.hstack((q_guess[:nd], mvi.lambda1, t_guess))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0,full_output=True, fprime=fprime)

    if ier != 1:
        # The root solver will fail for exactly simultaneous collisions since setting
        # t1 = t2 will return nan for the system's velocity. Can't handle this right now, 
        # but can at least check for it.
        if mvi.t1 == x[-1]:
            print 'Simultaneous impacts detected.'
        raise trep.ConvergenceError("Find impact failed to converge.")           

    mvi.q2 = np.append(x[:nd], kin_config(x[-1]))
    mvi.t2 = x[-1]
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lam1': mvi.lambda1, 'lam0': mvi.lambda0}


def plastic_impact(mvi):
    """ Finds the configuration after impact.
          when passed to this function, mvi should have following states
            0 - two states before impact
            1 - state just before impact
            2 - impact state
    """
    nd, m = mvi.nd, mvi.nc

    mvi.set_midpoint()
    tau = D4L2(mvi)

    mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    mvi.set_single_state(1, mvi.t2, mvi.q2, mvi.p2)
    mvi.t2 = mvi.t_end
    mvi.lambda0 = mvi.lambda1

    set_system_state(mvi, 1)
    Dh1T  = np.transpose(Dh(mvi))
    Dphi1 = np.transpose(Dphi(mvi))
    D2h1  = D2h(mvi) 

    def func(x):
        mvi.q2 = np.append(x[:nd], mvi.qk_end)
        mvi.lambda1 = x[nd:(nd+m)]
        lambda_c = x[-2]
        Ec = x[-1]        

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f4 = phi(mvi)

        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)-np.transpose(Dphi1)*lambda_c
        f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

        return np.hstack((f1, f2, f3, f4))

    def fprime(x):
        mvi.q2 = np.append(x[:nd], mvi.qk_end)
        mvi.lambda1 = x[nd:nd+m]
        lambda_c = x[-2]
        Ec = x[-1] 

        mvi.set_midpoint()
        Df11 = (D2D1L2(mvi) + D2fm2(mvi))
        Df31 = D2D3L2(mvi)

        set_system_state(mvi, 2)
        Df21 = Dh(mvi)
        Df41 = Dphi(mvi)
        Df32 = -D2h(mvi)

        Df1 = np.column_stack((Df11, -Dh1T, Dphi1, np.zeros((nd,1))))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc+2))))
        Df3 = np.hstack((Df31, Df32, (0.0, 1.0)))
        Df4 = np.hstack((Df41,  np.zeros(mvi.nc+2)))

        return np.vstack((Df1, Df2, Df3, Df4))

    x0 = np.hstack((mvi.q2[:nd], mvi.lambda1, np.dot(Dphi1, mvi.p1), 0.0))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Plastic impact update failed to converge.")  

    mvi.q2 = np.append(x[:nd],mvi.qk_end)
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1, 'lambdac': x[-2], 'Ec': x[-1]}

def elastic_impact(mvi, Ec):
    """ Finds the configuration after impact.
          when passed to this function, mvi should have following states
            0 - two states before impact
            1 - state just before impact
            2 - impact state
    """
    nd, m = mvi.nd, mvi.nc

    mvi.set_midpoint()
    tau = D4L2(mvi)

    mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    mvi.set_single_state(1, mvi.t2, mvi.q2, mvi.p2)
    mvi.t2 = mvi.t_end
    mvi.lambda0 = mvi.lambda1

    set_system_state(mvi, 1)
    Dh1T  = np.transpose(Dh(mvi))
    Dphi1 = np.transpose(Dphi(mvi))
    D2h1  = D2h(mvi) 

    def func(x):
        mvi.q2 = np.append(x[:nd], mvi.qk_end)
        mvi.lambda1 = x[nd:(nd+m)]
        lambda_c = x[-1]      

        set_system_state(mvi, 2)
        f2 = h(mvi)
 
        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)-np.transpose(Dphi1)*lambda_c
        f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

        return np.hstack((f1, f2, f3))

    def fprime(x):
        mvi.q2 = np.append(x[:nd], mvi.qk_end)
        mvi.lambda1 = x[nd:(nd+m)]
        lambda_c = x[-2]
        Ec = x[-1] 

        mvi.set_midpoint()
        Df11 = D2D1L2(mvi)+D2fm2(mvi)
        Df31 = D2D3L2(mvi)

        set_system_state(mvi, 2)
        Df21 = Dh(mvi)
        Df32 = -D2h(mvi)

        Df1 = np.column_stack((Df11, -Dh1T, Dphi1))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc+1))))
        Df3 = np.hstack((Df31, Df32, (0.0)))

        return np.vstack((Df1, Df2, Df3))

    x0 = np.hstack((mvi.q0[:nd], mvi.lambda1, np.dot(Dphi1, mvi.p1)))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Elastic impact update failed to converge.")  

    mvi.q2 = np.append(x[:nd],mvi.qk_end)
    mvi.lambda1 = x[nd:(nd+m)]
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

    # Perfectly inelastic.
    if mvi.cor == 0:
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

    # Perfectly elastic.
    elif mvi.cor == 1:
        elastic_impact(mvi, 0.0)
    
    # CoR between 0 and 1; need to solve twice.
    else:    
        next_state = plastic_impact(mvi)
        Ec_prime = mvi.cor*next_state['Ec']
        
        # At this point it might be nice to check if the energy lost is 
        # small enough to treat the impact as plastic.

        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.set_single_state(1, mvi.t0, mvi.q0, mvi.p0)
        mvi.lambda1 = mvi.lambda0

        elastic_impact(mvi, Ec_prime)              
