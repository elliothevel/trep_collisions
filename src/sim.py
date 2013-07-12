import trep

def simulate(cmvi, tf, dt, k2_func=None):
    """ Simulates the mvi to time tf with steps dt. This includes collisions. """
    
    # Extract initial conditions.
    q   = [cmvi.q2]
    t   = [cmvi.t2]
    lam = [cmvi.lambda1]

    # Simulation loop.
    while cmvi.t1 < tf:
        
        # Update dynamics. 

        # If there are kinematic variables, they should be calculated by 
        #   k2_func, which should accept the arguments (q, t) and return 
        #   a vector (tuple) of length mvi.nd.
        if k2_func is None:
            k2 = tuple()
        else:
            k2 = k2_func(cmvi.q1, cmvi.t1)
        
        try:            
            cmvi.step_(cmvi.t2+dt, k2=k2)
        except trep.ConvergenceError as e:
            print 'Convergence Error:', e.message
            print '\t\t   Ending simulaion at t=%f.' %cmvi.t1
            break    

        # Update the solution arrays.   
        q.append(cmvi.q2)
        t.append(cmvi.t2)
        lam.append(cmvi.lambda1)   

    # Return the time, configuration, and constraint force vectors.
    return (t, q, lam)
