import trep
import numpy as np
from detect import *
from impact import *
from release import *


class CollisionMVI(trep.MidpointVI):
    """ Collision midpoint variational integrator. """

    def __init__(self, system, surfaces, kin_func=None, tolerance=1e-10, 
                 impact_method='root', release_method='root', cor=0, energy_threshold=0.0):
        trep.MidpointVI.__init__(self, system, tolerance=tolerance)

        # Impact surfaces.
        self._surfaces = surfaces
        self._surface = None
                  
        # Method for solving release state.          
        assert release_method in ['root', 'interp', 'endpoint'], "Invalid release method"
        self._release_solve_method = release_method
        self._releases_off = False

        # Method for solving impact state.
        assert impact_method in ['root', 'interp', 'endpoint'], "Invalid impact method"
        self._impact_solve_method = impact_method

        # Coefficient of restitution for impacts and threshold for treating elastic 
        # impacts as plastic and for constraint releases.
        self.cor = cor
        self.energy_threshold = energy_threshold

        # Function for calculating kinematic configuration.
        self.kin_func = kin_func


    @property
    def surfaces(self):
        return self._surfaces

    @property
    def active_surfaces(self):
        return [surface for surface in self._surfaces if surface.active]

    @property
    def inactive_surfaces(self):
        return [surface for surface in self._surfaces if not surface.active]  

    @property
    def state0(self):
        return {'t': self.t0, 'q': self.q0, 'p': self.p0}

    @property
    def state1(self):
        return {'t': self.t1, 'q': self.q1, 'p': self.p1}

    @property
    def state2(self):
        return {'t': self.t2, 'q': self.q2, 'p': self.p2}

    @property
    def lambda1c(self):
        set_system_state(self, 1)
        return self.system.lambda_()

    @property
    def lambda2c(self):
        set_system_state(self, 2)
        return self.system.lambda_()

    def initialize_state(self, t1, q1, p1, lambda0=None, lambda1=None):
        """ Calls trep's built in initializer and also sets state 0. """
        self.initialize_from_state(t1, q1, p1, lambda1=lambda1)
        self.lambda0 = lambda0
        self.set_single_state(0, t1, q1, p1)
        self._dq1 = 0

    def set_single_state(self, num, t, q, p):
        """ Sets state 0, 1, or 2 of the integrator to a given (t,q,p). """
        if num == 0:
            self.t0, self.q0, self.p0 = t, q, p
        if num == 1:
            self.t1, self.q1, self.p1 = t, q, p
        if num == 2:
            self.t2, self.q2, self.p2 = t, q, p  

    def kin_config(self, t):
        if self.kin_func is None:
            return tuple()
        return self.kin_func(t)

    def cstep(self, t2, k2=tuple(), max_iterations=200, q2_hint=None, lambda1_hint=None):
        """ Steps the system to time t2. This satisfies the dynamics and resolves collisions. """

        # First have trep step the normal mvi.
        self.lambda0 = self.lambda1
        self.set_single_state(0, self.t1, self.q1, self.p1)
        if self.t1 != self.t2:
            self._dq1 = (self.q2-self.q1)/(self.t2-self.t1)
        self.step(t2, k2=self.kin_config(t2), q2_hint=q2_hint, lambda1_hint=lambda1_hint)

        # Then solve for a new state (t2, q2, p2) consistent with the impact surfaces.
        self.solve_collisions(self.t2)
        
    def solve_collisions(self, t2):
        """ Checks for impacts or releases and updates state 2 if necessary. """
        self.t_start, self.t_end = self.t1, t2
        self.current_dt = self.t_end-self.t_start

        # Loop until an acceptable state is found.
        while True:

            # Check if any events have occurred.
            impact_surfaces  = detect_impacts(self)
            release_surfaces = detect_releases(self)

            # Case 1: no events, return with the current state.
            if not any((release_surfaces, impact_surfaces)):
                return    

            # Case 2: impact but no releases, solve for impacts.
            elif (impact_surfaces and (not release_surfaces)): 
                surfaces, impact_state = solve_impact(self, impact_surfaces, self._impact_solve_method)
                impact_update(self, surfaces, impact_state) 
                                
            # Case 3: releases but no impacts, solve for releases.    
            elif (release_surfaces and (not impact_surfaces)):
                surfaces, release_state = solve_release(self, release_surfaces, self._release_solve_method)
                release_update(self, surfaces, release_state)  
 
            # Case 4: both impacts and releases, find first event.
            else:
                impact_surfaces, impact_state = solve_impact(self, impact_surfaces, self._impact_solve_method)
                release_surfaces, release_state = solve_release(self, release_surfaces, self._release_solve_method)

                # Times are very close to one another - transfer of contact.
                if abs(impact_state['t']-release_state['t'])<1e-3:
                    contact_transfer_update(self, impact_surfaces, release_surfaces, impact_state)

                # Impact occurred first.
                elif impact_state['t'] < release_state['t']:
                    impact_update(self, impact_surfaces, impact_state)

                # Release occurred first.
                else:
                    release_update(self, release_surfaces, release_state)

                    
    def prepare_to_visualize(self):
        """ Should be called before visualizing system. It will draw the constraint(s)
            even if there are no active surfaces at the end of the simulation. """
        self.system.hold_structure_changes()
        for surface in self.inactive_surfaces:
            surface.activate_constraint()
        self.system.resume_structure_changes()   
