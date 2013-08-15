""" I am removing support for elastic collisions completely because it simplifies handling
    multiple and instantaneous events. These are the remnants of all the elastic functions.
"""

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
   




def impact_update2(mvi, impact_surfaces, impact_state):
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
            print next_state['Ec']

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
      
