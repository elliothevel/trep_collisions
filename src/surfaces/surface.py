import trep


class Surface(trep.Constraint):
    def __init__(self, system, name, tolerance=1e-10):
        system.hold_structure_changes()
        super(Surface, self).__init__(system, name=name, tolerance=tolerance) 
        self.system._constraints = self.system._constraints[:-1]
        self.system.resume_structure_changes()
        self.active = False
        
    def deactivate_constraint(self):
        self.system._constraints = tuple([con for con in self.system.constraints if (con != self)]) # TUPLE!!
        self.active = False 

    def activate_constraint(self):
        self.system._add_constraint(self)
        self.active = True

    def h(self):
        return self.phi()

    def h_dq(self, q1):
        return self.phi_dq(q1)

    def h_dqdq(self, q1, q2):
        return self.phi_dqdq(q1, q2)

    def h_dqdqdq(self, q1, q2, q3):
        return self.phi_dqdqdq(q1, q2, q3)


def global_surfaces(surface, system, impact_frames=[], tolerance=1e-10, kwargs={}):
    if not impact_frames: 
        impact_frames = system.masses
    return [surface(system, frame=frame1, tolerance=tolerance, **kwargs) for frame1 in impact_frames]


