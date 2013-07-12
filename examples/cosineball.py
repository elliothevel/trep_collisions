import trep
from trep import ty, tz, rx
from trep.visual import *
import numpy as np

import trep_collisions as tc

# Create the system.
ball = trep.System()
ball.import_frames([ty('y'), [tz('z', mass=1.0)]])
trep.potentials.Gravity(ball)

# Create the surface.
surface = tc.surfaces.Cosine(ball, 1, 1)

# Initialize the integrator.
mvi = tc.CollisionMVI(ball, surface)

# Give the initial condition.
ball.get_config('z').q = 7.0
ball.get_config('y').q = 0.7
mvi.initialize_state(0.0, ball.q, tc.util.D2L2(mvi))

# Simulate.
tf = 5.0
dt = 0.01
(t, q, lam) = tc.simulate(mvi, tf, dt)

# Visualize.
mvi.prepare_to_visualize()
visualize_3d([VisualItem3D(ball, t, q)])
