import trep
import numpy as np
import scipy
import scipy.optimize
from detect import *
from util import *
from release import *

TIME_TOL = 1e-7


def linear_interpolate_impact(mvi):
    """ Finds a time and configuration when phi=0 based on linear interpolation.
        this is used to get a guess to seed the root-solver. """

    set_system_state(mvi, 1)
    phi0 = phi(mvi)
    set_system_state(mvi, 2)
    phi1 = phi(mvi)

    t_i = mvi.t1 - (phi0/(phi1-phi0))*(mvi.t2-mvi.t1)
    frac = (t_i-mvi.t1)/(mvi.t2-mvi.t1)
    q_i = frac*(mvi.q2-mvi.q1)+mvi.q1

    return t_i, q_i


def solve_impact_time(mvi):
    """ Finds the time and configuration at impact.
        when passed to this function, mvi should have following states:
            0 - two states before impact
            1 - state just before impact
            2 - invalid configuration to be corrected

        solve_impact_time modifies only state 2, lambda1    
    """

    nd, m = mvi.nd, mvi.nc

    set_system_state(mvi, 1)
    Dh1T = np.transpose(Dh(mvi))

    # Here, check if the impact occurs exactly at t1. If phi is within tolerance, set the
    # impact state at 1. The root-solver can't find this because it will cause a divide by zero 
    # error.
    if abs(phi(mvi))<mvi._surface.tolerance:
        mvi.q2 = mvi.q1
        mvi.t2 = mvi.t1
        mvi.lambda1 = mvi.lambda0   
        mvi.p2 = mvi.p1
        return {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lambda1': mvi.lambda1}

    def func(x):  
        mvi.t2 = x[-1]
        mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))
        mvi.lambda1 = x[nd:(nd+m)]

        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f3 = phi(mvi)

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
        Df31 = Dphi(mvi)
        Df23 = D2h(mvi).reshape((mvi.nc, 1))

        Df1 = np.hstack((Df11, -Dh1T, Df13))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc)), Df23))
        Df3 = np.hstack((Df31, np.zeros(mvi.nc+1)))

        return np.vstack((Df1, Df2, Df3))

    
    t_guess, q_guess = linear_interpolate_impact(mvi)
    x0 = np.hstack((q_guess[:nd], mvi.lambda1, t_guess))
    try:
        (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)
    except ZeroDivisionError: 
        raise Exception("Divide by zero in impact time solver.")

    if ier != 1:
        raise trep.ConvergenceError("Find impact failed to converge.") 

    mvi.q2 = np.append(x[:nd], mvi.kin_config(x[-1]))
    mvi.t2 = x[-1]
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)
    return {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}


def plastic_impact(mvi, impact_surfaces):
    """ Finds the configuration after impact.
          when passed to this function, mvi should have following states
            0 - two states before impact
            1 - state just before impact
            2 - impact state

        plastic_impact updates so that state2 is post-impact, state1 is impact, state0 is before impact.       
    """
    p = len(impact_surfaces)
    nd, m = mvi.nd, mvi.nc

    # Special case: impact occurs at t1. Step integrator back.
    if mvi.t2 == mvi.t1:
        mvi.t1 = mvi.t0
        mvi.q1 = mvi.q0
        mvi.p1 = mvi.p0

    mvi.set_midpoint()
    tau = D4L2(mvi)

    mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    mvi.set_single_state(1, mvi.t2, mvi.q2, mvi.p2)
    mvi.lambda0 = mvi.lambda1

    # It is always possible that an impact happens very close to the 
    # time step t2. If this is the case, then we need to take a whole new
    # step to ensure the impact map is computed. This step is assumed to 
    # be the same size as the previous one. Not ideal, but definitely needed.
    if (mvi.t2-mvi.t_end)<TIME_TOL:
        mvi.t_end = mvi.t_start+2*mvi.current_dt 
    mvi.t2 = mvi.t_end

    set_system_state(mvi, 1)
    Dh1T  = np.transpose(Dh(mvi))
    D2h1  = D2h(mvi) 
    Dphi1 = Dphi(mvi, surfaces=impact_surfaces)

    def func(x):
        mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
        mvi.lambda1 = x[nd:(nd+m)]
        lambda_c = x[(nd+m):-1]
        Ec = x[-1]        

        set_system_state(mvi, 2)
        f2 = h(mvi)
        f4 = phi(mvi, surfaces=impact_surfaces)
            
        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)-np.dot(np.transpose(Dphi1),lambda_c)
        f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

        return np.hstack((f1, f2, f3, f4))


    def fprime(x):
        mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
        mvi.lambda1 = x[nd:nd+m]
        lambda_c = x[(nd+m):-1]
        Ec = x[-1] 

        mvi.set_midpoint()
        Df11 = D2D1L2(mvi) + D2fm2(mvi)
        Df31 = D2D3L2(mvi)

        set_system_state(mvi, 2)
        Df21 = Dh(mvi)
        Df41 = Dphi(mvi, surfaces=impact_surfaces)
        Df32 = -D2h(mvi)

        Df1 = np.column_stack((Df11, -Dh1T, -np.transpose(Dphi1), np.zeros((nd,1))))
        Df2 = np.hstack((Df21, np.zeros((mvi.nc, mvi.nc+p+1))))
        Df3 = np.hstack((Df31, Df32, np.zeros(p), (1.0,)))
        Df4 = np.hstack((Df41,  np.zeros((p,mvi.nc+p+1))))

        return np.vstack((Df1, Df2, Df3, Df4))


    x0 = np.hstack((mvi.q2[:nd], mvi.lambda1, np.dot(Dphi1, mvi.p1), np.dot(D2h1, mvi.lambda1)-tau))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if (ier != 1) and scipy.linalg.norm(func(x))>mvi.tolerance:
        print mesg  
        raise trep.ConvergenceError("Plastic impact update failed to converge.")  

    mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1, 'lambdac': x[(nd+m):-1], 'Ec': x[-1]}


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
    mvi.lambda0 = mvi.lambda1
    if (mvi.t2-mvi.t_end)<TIME_TOL:
        mvi.t2 = mvi.t_end+(mvi.t_end-mvi.t_start)
    else:    
        mvi.t2 = mvi.t_end

    set_system_state(mvi, 1)
    Dh1T  = np.transpose(Dh(mvi))
    Dphi1 = np.transpose(Dphi(mvi))
    D2h1  = D2h(mvi) 

    def func(x):
        mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
        mvi.lambda1 = x[nd:(nd+m)]
        lambda_c = x[-1]      

        set_system_state(mvi, 2)
        f2 = h(mvi)
 
        mvi.set_midpoint()
        f1 = mvi.p1+D1L2(mvi)+fm2(mvi)-np.dot(Dh1T,mvi.lambda1)-np.transpose(Dphi1)*lambda_c
        f3 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

        return np.hstack((f1, f2, f3))

    def fprime(x):
        mvi.q2 = np.append(x[:nd], mvi.kin_config(mvi.t2))
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

    x0 = np.hstack((mvi.q1[:nd], mvi.lambda1, np.dot(Dphi1, mvi.p1)))
    (x, infodict, ier, mesg) = scipy.optimize.fsolve(func, x0, full_output=True, fprime=fprime)

    if ier != 1:
        print mesg
        raise trep.ConvergenceError("Elastic impact update failed to converge.")  

    mvi.q2 = np.append(x[:nd],mvi.kin_config(mvi.t2))
    mvi.lambda1 = x[nd:(nd+m)]
    mvi.set_midpoint()    
    mvi.p2 = D2L2(mvi)

    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1, 'lambdac': x[-2], 'Ec': x[-1]}



def simultaneous_elastic_impact(mvi, Ec, surfaces):
    """ Finds the configuration after impact.
          when passed to this function, mvi should have following states
            0 - two states before impact
            1 - state just before impact
            2 - impact state

    Here the unknown variables are dq and lambda for contact surfaces.
            
    """
    print 'energy', Ec
    nd, m = mvi.nd, mvi.nc
    p = len(surfaces)

    mvi.set_midpoint()
    tau = D4L2(mvi)

    #mvi.set_single_state(0, mvi.t1, mvi.q1, mvi.p1)
    #mvi.set_single_state(1, mvi.t2, mvi.q2, mvi.p2)

    set_system_state(mvi, 2)
    Dh1T  = np.transpose(Dh(mvi))
    Dphi1T = np.transpose(Dphi(mvi, surfaces=surfaces))
    D2h1  = D2h(mvi) 

    def func(x):         
        mvi.system.set_state(q=mvi.q2, dq=x[:nd], t=mvi.t2)
        f1 = mvi.p1+D1L2_lim(mvi)-np.dot(Dphi1T,x[(nd+m):])
        f2 = tau+D3L2(mvi)-np.dot(D2h1, mvi.lambda1)+Ec

        #print 'f', np.linalg.norm(np.hstack((f1, f2)) )
        return np.linalg.norm(np.hstack((f1, f2)))

    x0 = np.hstack((mvi.system.dq[:nd], np.dot(np.transpose(Dphi1T), mvi.p2)))
    x, fx,its, imode, smode = scipy.optimize.fmin_slsqp(func, x0, full_output=True, acc=1e-8)


    if imode != 0:
        print smode
        raise trep.ConvergenceError("Simultaneous elastic impact update failed to converge.")  
    print 'solved', x0[:nd], x[:nd]

    mvi.system.set_state(q=mvi.q2, dq=x[:nd], t=mvi.t2)    
    mvi.p2 = D2L2_lim(mvi)

    # At this point we have the configuration and momentum at the point of impact.
    # We finish by simulating to t2.
    mvi.lambda0 = mvi.lambda1
    if abs(mvi.t_end-mvi.t2)<TIME_TOL:
        t2 = mvi.t_start+mvi.current_dt
    else:
        t2 = mvi.t_end

    mvi.step(t2, mvi.kin_config(t2))
    print detect_releases(mvi)
    print detect_impacts(mvi)
    
    return {'t': mvi.t2, 'q': mvi.q2, 'lambda1': mvi.lambda1}
    
 
def solve_impact(mvi, impact_surfaces, method):
    """ When given one or more impacts, finds the first impact to occur and the state
        at that point. Returns the frame and state.
        
        solve impact returns mvi with the same states as it was passed with.    
    """
        
    # Endpoint method: all impacts occur at state 1. 
    if method == 'endpoint':
        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.lambda1 = mvi.lambda0
        state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}
        return impact_surfaces, state

    # Interp method: linearly interpolate impact state.
    elif method == 'interp':
        times = []
        for surf in impact_surfaces:
            mvi._surface = surf
            times.append(linear_interpolate_impact(mvi)[0])
        ti = np.mean(times)
        mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
        mvi.step(ti, k2=mvi.kin_config(ti))
        state = {'t': mvi.t2, 'q': mvi.q2, 'p': mvi.p2, 'lambda1': mvi.lambda1}
        return impact_surfaces, state

    # Root method: root-solve the time of impact. Solve separately and then
    # check for multiple impacts.
    elif method == 'root':

        q2, p2, t2, lam1 = mvi.q2, mvi.p2, mvi.t2, mvi.lambda1

        # First solve for each impact separately, saving the states.
        impacts = {}
        for surface in impact_surfaces:
            mvi._surface = surface
            s = solve_impact_time(mvi)
            impacts.update({surface: s})
            mvi.q2, mvi.p2, mvi.t2, mvi.lambda1 = q2, p2, t2, lam1

        # Find the first surface to impact assuming distinct impacts.
        impact_surface = impacts.keys()[np.argmin([impacts[k]['t'] for k in impacts.keys()])]
        impact_state = impacts[impact_surface]

        # Finally, check for multiple impacts. Multiple impacts are said to have occurred
        # if the impact configuration lies on more that one surface.
        surfaces = [impact_surface]
        mvi.system.set_state(q=impact_state['q'], t=impact_state['t'])
        for surface in impact_surfaces:
            if surface != impact_surface:
                if abs(surface.phi())<surface.tolerance:
                    surfaces.append(surface)         

        return surfaces, impact_state   


def impact_update(mvi, impact_surfaces, impact_state):
    """ Given a state at the time of impact, this function solves the impact map and 
        adds the new constraint to the system.
        
        when passed to impact_update, state2 is invalid    
    """  
     
    # Set system config to impact config. Solve impact update.
    mvi.set_single_state(2, impact_state['t'], impact_state['q'], impact_state['p'])
    mvi.lambda1 = impact_state['lambda1']

    # Here we need to make sure no constraints are released in the process of the impact update.
    # If there are, we need to deal with them here because solving them later will use a regular
    # integrator step instead of the impact map.
    while True:
        treat_as_plastic = False

        # Perfectly plastic impact.
        if mvi.cor == 0:
            next_state = plastic_impact(mvi, impact_surfaces) 

        # Perfectly elastic impact.
        elif mvi.cor == 1:
            if len(impact_surfaces)==1:
                elastic_impact(mvi, 0.0)
            else:
                print 'Warning: simultaneous elastic impact.'
                simultaneous_elastic_impact(mvi, 0.0, impact_surfaces)
    
        # Impact with cor in (0,1).
        else:    
            next_state = plastic_impact(mvi, impact_surfaces)
            Ec_prime = (1-mvi.cor)*next_state['Ec']
            print next_state
        
            mvi.set_single_state(2, mvi.t1, mvi.q1, mvi.p1)
            mvi.set_single_state(1, mvi.t0, mvi.q0, mvi.p0)
            mvi.lambda1 = mvi.lambda0

            # At this point we check if the energy lost is 
            # small enough to treat the impact as plastic.
            if Ec_prime < mvi.energy_threshold:
                next_state = plastic_impact(mvi, impact_surfaces)
                treat_as_plastic = True
            else: 
                if len(impact_surfaces)==1:
                    elastic_impact(mvi, Ec_prime)  
                else:
                    print 'Warning: simultaneous elastic impact.'
                    simultaneous_elastic_impact(mvi, Ec_prime, impact_surfaces)

        if mvi._release_solve_method != 'endpoint':
            break

        else:
            release_check = detect_releases(mvi)

            if not release_check:
                break

            else:
                release_update(mvi, release_check, {'t': mvi.t1, 'q': mvi.q1, 'p': mvi.p1, 'lambda0': mvi.lambda0}, stepafter=False)


     # Perfectly inelastic, add constraint to system.
    if (mvi.cor == 0) or treat_as_plastic:
        state = save_state(mvi)  

        mvi.system.hold_structure_changes()
        for surface in impact_surfaces:
            surface.activate_constraint()
        mvi.system.resume_structure_changes()  
        #mvi._structure_updated()

        # Adding constraints can sometimes clear the inteagrator's states. Here, we restore the
        #   current state. 
        restore_state(mvi, state)
        mvi.lambda0 = np.append(mvi.lambda0, [0.0]*len(impact_surfaces)) 
        mvi.lambda1 = np.append(next_state['lambda1'],  next_state['lambdac'])
       
