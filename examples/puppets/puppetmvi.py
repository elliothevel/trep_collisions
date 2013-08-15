import trep_collisions as tc
import time


class PuppetMVI(tc.CollisionMVI):
    """ A type of collision integrator that handles string constraints
        in a special way.
    """

    def __init__(self, system, collision_surfaces=[], puppet_strings=[], **kwargs):
        surfaces = collision_surfaces + puppet_strings
        tc.CollisionMVI.__init__(self, system, surfaces, **kwargs)
        self._strings = puppet_strings
        self._start = time.time()

    @property
    def active_strings(self):
        return [string for string in self._strings if string.active]

    @property
    def inactive_strings(self):
        return [string for string in self._strings if not string.active] 

    def release_all_strings(self):
        tc.release.remove_constraints(self, self.active_strings)

    def add_all_strings(self):
        tc.impact.add_constraints(self, self.inactive_strings)

    def check_string_lengths(self):
        return [s.h() for s in self.inactive_strings]

    def release_all_contacts(self):
        contacts = [surf for surf in self.active_surfaces if (surf not in self._strings)]
        tc.release.remove_constraints(self, contacts)

    def cstep(self, *args, **kwargs):
        start = time.time()
        super(PuppetMVI, self).cstep(*args, **kwargs)
        print 't=%f, steptime=%f, cumtime=%f' %(self.t2, time.time()-start, time.time()-self._start)

    def impact_update(self, impact_surfaces, impact_state):
        """ This method overrides the one from the normal CMVI. The difference is that the 
            string constraints are by default deactivated for the impact updates.
        """   
        self.release_all_strings()

        # Set system config to impact config.
        self.set_single_state(2, impact_state['t'], impact_state['q'], impact_state['p'])
        self.tau2 = impact_state['tau']

        # Apply the impact map.
        next_state = tc.impact.plastic_impact(self, impact_surfaces) 

        # Add the constraints to the system.
        tc.impact.add_constraints(self, impact_surfaces)

        # Add strings back.
        add_strings = []
        for s in self._strings:
            if s not in impact_surfaces:
                if s.h()<0:
                    add_strings.append(s)
        tc.impact.add_constraints(self, add_strings)      
    
