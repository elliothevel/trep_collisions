import trep_collisions as tc
from strings import *


class PuppetMVI(tc.CollisionMVI):
    """ This is a midpoint variational integrator + collision solver
        that also handles puppet strings. 
    """

    def __init__(self, system, surface, impact_frames=None, puppet_strings=[], tolerance=1e-10, cor=0.0):
        tc.CollisionMVI.__init__(self, system, surface, impact_frames, 
                                 tolerance=tolerance, root_solve=False, cor=0.0)

        self.puppet_strings = puppet_strings
            
        if impact_frames is None:
            self.impact_frames = system.masses
        else:
            self.impact_frames = impact_frames  
              
    def step_(self, t2, k2=tuple(), max_iterations=200, q2_hint=None,
              lambda1_hint=None, verbose=False):
        """ Steps the system from t1 --> t2. """
        # First have trep step the normal mvi.
        self.lambda0 = self.lambda1
        self.t0 = self.t1
        self.q0 = self.q1
        self.p0 = self.p1
        self.step(t2, k2=k2, q2_hint=q2_hint, lambda1_hint=lambda1_hint)

        # Then solve for a new state (q2, p2) consistent with the impact surface.
        self.solve_collisions(self.t2)
    

    def solve_collisions(self, t2):
        """ Checks for impacts or releases and updates state 2 if necessary. """
        self.t_end = t2  
        self.qk_end = self.q2[self.nd:]  

        # Loop until an acceptable state is found.
        while True:

            # Check if any events have occurred.
            impact_frames   = detect_impacts(self)
            release_frames  = detect_releases(self)
            catch_strings   = detect_string_cathces(self)
            release_strings = detect_string_releases(self)

            # Case 1: no events, return with the current state.
            if  ((not release_frames) and (not impact_frames)
                 and (not catch_strings) and (not release_strings)):
                return    

            # Case 2: One or more releases.    

            # For the PuppetMVI, release root-solving is turned off. This means all releases will
            # occur at t1. Thus, if we have any releases, handle them first.
            if (release_frames or release_strings):
                release_update_with_strings(self, release_frames, release_strings, self.t1, self.q1, self.p1, self.lambda0)


            # Case 3: No releases, string catches or impacts.
            elif (impact_frames and (not catch_strings)):
                frame, impact = solve_impact(self, impact_frames)
                impact_update(self, frame, impact['t'], impact['q'], impact['p'], impact['lam0'], impact['lam1'])
                
            elif (catch_strings and (not impact_frames)):
                catch_time = get_catch_time(self, catch_strings)
                solve_string_catch(self, catch_strings, catch_time)
                
            else:
                frame, impact = solve_impact(self, impact_frames)
                catch_time = get_catch_time(self, catch_strings)

                if catch_time <= impact['t']:
                    solve_string_catch(self, catch_strings, catch_time)

                else:
                    impact_update(self, frame, impact['t'], impact['q'], impact['p'], impact['lam0'], impact['lam1'])    

