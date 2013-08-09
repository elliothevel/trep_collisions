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

# Create the surface.
surface = tc.surfaces.global_surfaces(tc.surfaces.Ground, ball)

# Initialize the integrator. Here we choose a coefficient of restitution
# between 0 and 1 to get elastic impacts with energy loss. The energy threshold
# specifies that if a collision results in an energy loss less than some value,
# the impact will be treated as plastic. This will result in the ball resting on
# the ground after several bounces.
mvi = tc.CollisionMVI(ball, surface, cor=0.5, energy_threshold=0.03)

# Give the initial condition.
ball.get_config('z').q = 5.0
ball.get_config('y').dq = 0.0
mvi.initialize_state(0.0, ball.q, tc.util.D2L2(mvi))

# Simulate.
tf = 10.0
dt = 0.01
(t, q, lam) = tc.simulate(mvi, tf, dt)

# Plot the height of the ball vs. time.
#h = [q1[ball.get_config('z').index] for q1 in q]
#plt.plot(t, h)
#plt.show()

# Visualize.
mvi.prepare_to_visualize()
visualize_3d([VisualItem3D(ball, t, q)])
