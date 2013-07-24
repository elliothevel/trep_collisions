"""
In this example, we have two pendula that can contact one another (basically newton's cradle
with only two masses). To solve this we show how a new collision surface can be defined and used.
"""

import trep
from trep import ty, tz, rx
from trep.visual import *
import trep_collisions as tc

# Define the system. We include gravity and damping.
frames = [ty(0.165), [rx('theta-r'), [tz(-1.0, mass=1.0)]],
          ty(-0.165), [rx('theta-l'), [tz(-1.0, mass=1.0)]]]
s = trep.System()
s.import_frames(frames)
trep.potentials.Gravity(s)
trep.forces.Damping(s, 0.25)

# Here we define a new collision surface so that the masses 
# contact one another.
class NewSurface:
    def __init__(self, system, frame1, frame2, dist):
        self.system = system
        self.frame1 = frame1
        self.frame2 = frame2
        self.tape_measure = trep.TapeMeasure(system, [frame1, frame2])
        assert dist != 0, "distance can't be 0"
        self.dist = dist

    def set_frame(self, frame):
        pass

    def phi(self):
        return self.tape_measure.length()-self.dist

    def phi_dq(self, q_i):
        return self.tape_measure.length_dq(q_i)

    def add_constraint_to_system(self, frame):
        trep.constraints.Distance(self.system, self.frame1, self.frame2, 
                                  self.dist, name='constraint_'+frame.name)

# Create an instance of the new surface and the integrator. Give initial
# conditions. We set the coefficient of restitution to be 0.5.
surf = NewSurface(s, s.masses[0], s.masses[1], 0.3)
cmvi = tc.CollisionMVI(s, surf, impact_frames=[s.masses[0]], cor=0.5)
cmvi.initialize_state(0.0, [1.0,-1.0], [0,0])

# Simulate the dynamics.
tf = 7.0
dt = 0.01
(t, q, lam) = tc.simulate(cmvi, tf, dt)

# Visulize the system.
cmvi.prepare_to_visualize()
visualize_3d([VisualItem3D(s, t, q)])
