"""
In this example, we have two pendula that can contact one another (basically newton's cradle
with only two masses).
"""

import trep
from trep import ty, tz, rx
from trep.visual import *
import numpy as np
import trep_collisions as tc

# Define the system. We include gravity and damping.
frames = [ty(0.165), [rx('theta-r'), [tz(-1.0, mass=1.0)]],
          ty(-0.165), [rx('theta-l'), [tz(-1.0, mass=1.0)]]]
s = trep.System()
s.import_frames(frames)
trep.potentials.Gravity(s)
trep.forces.Damping(s, 0.0)

# Create an instance of the new surface and the integrator. Give initial
# conditions. We set the coefficient of restitution to be 0.5.
surfaces = [tc.surfaces.Distance(s, s.masses[0], s.masses[1], 2*0.165, invalid='short')]
#surfaces = []
cmvi = tc.CollisionMVI(s, surfaces, release_method='interp', cor=0.0)
cmvi.initialize_state(0.0, [1.0,-1.5], [0,0])

# Simulate the dynamics.
tf = 2.0
dt = 0.01
(t, q, lam) = tc.simulate(cmvi, tf, dt)

# Visulize the system.
cmvi.prepare_to_visualize()
visualize_3d([VisualItem3D(s, t, q)])
