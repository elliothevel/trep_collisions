import trep
from trep import ty, tz, rx
from trep.visual import *
import numpy as np
import matplotlib.pyplot as plt

import trep_collisions as tc

# Create the system.
ball = trep.System()
ball.import_frames([ty('y'), [tz('z', mass=1.0)]])
trep.potentials.Gravity(ball)
#trep.forces.Damping(ball, 1.0)

# Create the surface.
surface = tc.surfaces.Ground(ball)

# Initialize the integrator.
mvi = tc.CollisionMVI(ball, surface, cor=0.5)

# Give the initial condition.
ball.get_config('z').q = 5.0
ball.get_config('y').dq = 0.0
mvi.initialize_state(0.0, ball.q, tc.util.D2L2(mvi))

# Simulate.
tf = 5.0
dt = 0.01
(t, q, lam) = tc.simulate(mvi, tf, dt)

# Plot the height of the ball vs. time.
#h = [q1[ball.get_config('z').index] for q1 in q]
#plt.plot(t, h)
#plt.show()

# Visualize.
mvi.prepare_to_visualize()
visualize_3d([VisualItem3D(ball, t, q)])
