import trep
import numpy as np
from detect import *
from impact import *
from release import *

class CollisionMVI(trep.MidpointVI):
    """ This is a midpoint variational integrator + collision solver. """

    def __init__(self, system, surface, impact_frames=None, tolerance=1e-10, 
                       root_solve=False, cor=0, energy_threshold=0.0):
        trep.MidpointVI.__init__(self, system, tolerance=tolerance)
        self.constrained_frames = []
        self.surface = surface
            
        if impact_frames is None:
            self.impact_frames = system.masses
        else:
            self.impact_frames = impact_frames    
        self.initialize_names()
        
        # The root-solver is dynamically accurate but frequently fails to converge.
        #   By default it is turned off, meaning the time of release is mvi.t1 when lambda
        #   changes sign.
        self.root_solve_release = root_solve

        # Coefficient of restitution for impact.
        self.cor = cor

        # Threshold for treating elastic impacts as plastic.
        self.energy_threshold = energy_threshold


    def initialize_names(self):
        """ Creates unique names for all frames that don't already have one. """ 
        i=0
        for frame in self.impact_frames:
            if frame.name is None:
                frame.name = 'frame_%d' %i
                i+=1
                

    def initialize_state(self, t1, q1, p1, lambda0=None, lambda1=None):
        """ Call's trep's built in initializer and also sets state 0. """
        self.initialize_from_state(t1, q1, p1, lambda1=lambda1)
        self.lambda0 = lambda0
        self.t0 = t1
        self.q0 = q1
        self.p0 = p1


    def set_single_state(self, num, t, q, p):
        """ Sets state 0, 1, or 2 of the integrator to a given (t,q,p). """
        if num == 0:
            self.t0, self.q0, self.p0 = t, q, p
        if num == 1:
            self.t1, self.q1, self.p1 = t, q, p
        if num == 2:
            self.t2, self.q2, self.p2 = t, q, p  


    def step_(self, t2, k2=tuple(), max_iterations=200, q2_hint=None, lambda1_hint=None):
        """ Steps the system to time t2. This satisfies the dynamics and resolves collisions. """

        # First have trep step the normal mvi.
        self.lambda0 = self.lambda1
        self.t0 = self.t1
        self.q0 = self.q1
        self.p0 = self.p1
        self.step(t2, k2=k2, q2_hint=q2_hint, lambda1_hint=lambda1_hint)

        # Then solve for a new state (t2, q2, p2) consistent with the impact surface.
        self.solve_collisions(self.t2)
        

    def solve_collisions(self, t2):
        """ Checks for impacts or releases and updates state 2 if necessary. """
        self.t_end = t2  
        self.qk_end = self.q2[self.nd:]  

        # Loop until an acceptable state is found.
        while True:

            # Check if any events have occurred.
            impact_frames  = detect_impacts(self)
            release_frames = detect_releases(self)

            # Case 1: no events, return with the current state.
            if  ((not release_frames) and (not impact_frames)):
                return    

            # Case 2: impact but no releases, solve for impacts.
            elif (impact_frames and (not release_frames)): 
                frame, impact = solve_impact(self, impact_frames)
                impact_update(self, frame, impact['t'], impact['q'], impact['p'], impact['lam0'], impact['lam1']) 
                                
            # Case 3: releases but no impacts, solve for releases.    
            elif (release_frames and (not impact_frames)):
                frame, release = solve_release(self, release_frames, root_solve=self.root_solve_release)
                release_update(self, frame, release['t'], release['q'], release['p'], release['lam0'])  
 
            # Case 4: both impacts and releases, find first event.
            else:
                if not self.root_solve_release:
                    release_update(self, release_frames, self.t1, self.q1, self.p1, self.lambda0) 

                else:
                    iframe, impact = solve_impact(self, impact_frames)
                    rframe, release = solve_release(self, release_frames)

                    # Impact occurred first.
                    if impact['t'] < release['t']:
                        impact_update(self, iframe, impact['t'], impact['q'], impact['p'],
                                                    impact['lam0'], impact['lam1'])

                    # Release occurred first.
                    else:
                        release_update(self, rframe, release['t'], release['q'], release['p'], release['lam0'])
                    

    def prepare_to_visualize(self):
        """ Should be called before visualizing system. It will draw the constraint(s)
            even if there are no constrained frames at the end of the simulation. """
        if len(self.constrained_frames) == 0:
            self.surface.add_constraint_to_system(self.system.world_frame)  
