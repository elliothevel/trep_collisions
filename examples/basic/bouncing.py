import trep
from trep import ty, tz, rx
from trep.visual import *
import numpy as np
import matplotlib.pyplot as plt

import trep_collisions as tc


class Triangle(trep.System):
    def __init__(self, m, a):
        super(Triangle, self).__init__()

        self.import_frames(self.make_frames(m, a))
        #self.corners =  [self.get_frame('p%d'%i) for i in range(1,4)]
        trep.potentials.Gravity(self)

    def make_frames(self, m, a):
        return [ty('y'), [tz('z', mass=0.2*m, name='center'), 
                [rx('theta', mass=(0, m*a**2/6, 0, 0)), 
                 [rx(0.0), [tz(a, mass=0.2*m)],
                  rx(np.radians(120)), [tz(a, mass=0.2*m)],
                  rx(np.radians(240)), [tz(a, mass=0.2*m)]]]]]


# Create the system.
#ball = trep.System()
#ball.import_frames([ty('y1'), [tz('z1', mass=1.0, name='1'),
#                    [rx(1.0), [tz(-1.0, mass=1.0, name='2')]]]]),
#trep.potentials.Gravity(ball)
#trep.constraints.Distance(ball, '1', '2', 1.0)
#trep.constraints.Distance(ball, '2', '3', 1.0)

#trep.potentials.LinearSpring(ball, '1', '2', 10.0)
ball = Triangle(1.0, .50)

# Create the surface.
surface = tc.surfaces.global_surfaces(tc.surfaces.Ground, ball, impact_frames=ball.masses)
#surface = []

# Initialize the integrator. Here we choose a coefficient of restitution
# between 0 and 1 to get elastic impacts with energy loss. The energy threshold
# specifies that if a collision results in an energy loss less than some value,
# the impact will be treated as plastic. This will result in the ball resting on
# the ground after several bounces.
mvi = tc.CollisionMVI(ball, surface, cor=0.0, energy_threshold=0.03, release_method='interp')
#mvi._releases_off = True

# Give the initial condition.
ball.get_config('z').q = 3.0
ball.get_config('y').q = 0.0
ball.get_config('theta').q = 0.35
mvi.initialize_from_configs(0.0, (0,2.0,0.3), 0.01, (0.01,2.1,0.4))
#ball.satisfy_constraints()
mvi.initialize_state(0.0, mvi.q2, mvi.p2)

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
#visualize_3d([VisualItem3D(ball, [0, 1], [ball.q, ball.q])])
