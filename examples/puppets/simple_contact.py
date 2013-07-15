import trep
import trep_collisions as tc
import numpy as np
from puppet import *
from puppetmvi import *

# Create the puppet system, specifying which sets of strings should be used.
strings = {'shoulders': True, 'knees': False, 'hands': True}
puppet = Puppet(strings)

# Define the ground as the contact surface.
surface = tc.surfaces.Ground(puppet, dim=3)

# Initiate the integrator.
cmvi = tc.CollisionMVI(puppet, surface, impact_frames=puppet.feet_frames)
cmvi.initialize_state(0.0, puppet.q, np.zeros(cmvi.nd))

# Define a function that will move the puppet up and down. This is keeping the
# strings a constant length and moving the overhead frame up and down.
amp, freq = -1.0, 0.25
h0 = puppet.get_config('Frame Z').q
def h_func(q, t):
    h = h0 + amp*np.sin(freq*2*np.pi*t)
    return (h, )

# Simulate the dynamics.
tf = 2.5
dt = 0.01
(t, q, lam) = tc.simulate(cmvi, tf, dt, h_func)

# Visualize the results.
pos = [12.8121004 , -10.37581321, 7.70056414] 
ang = [2.42159265, 0.06 , 0.]
cmvi.prepare_to_visualize()
trep.visual.visualize_3d([PuppetVisual(puppet, t, q)], camera_pos=pos, camera_ang=ang)    
