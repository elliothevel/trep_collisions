"""
This demo simulates a puppet contacting the ground. The strings are modeled as unilateral
distance constraints with lengths controlled by kinematic configuration variables.
"""

import trep
import trep_collisions as tc
import numpy as np
from puppet import *
from puppetmvi import *

# Create the puppet system, specifying which sets of strings should be used.
strings = {'shoulders': True, 'knees': False, 'hands': True}
puppet = Puppet(strings)

# Define the ground as the contact surface.
collision_surfaces = tc.surfaces.global_surfaces(tc.surfaces.Ground, puppet, kwargs={'dim': 3}, impact_frames=puppet.feet_frames[:2])
puppet_strings = puppet.strings

# Define a function to move the puppet truck up and down.
rs0 = puppet.get_config('RightShoulderString').q
ls0 = puppet.get_config('LeftShoulderString').q
rh0 = puppet.get_config('RightHandString').q
lh0 = puppet.get_config('LeftHandString').q

rs_index = puppet.get_config('RightShoulderString').k_index
ls_index = puppet.get_config('LeftShoulderString').k_index
rh_index = puppet.get_config('RightHandString').k_index
lh_index = puppet.get_config('LeftHandString').k_index

amp, freq = 1.5, 0.2
k2 = np.zeros(len(puppet_strings))
def h_func(t):
    k2[rs_index] = rs0 + amp*np.sin(freq*2*np.pi*t)
    k2[ls_index] = ls0 + amp*np.sin(freq*2*np.pi*t)
    k2[rh_index] = rh0 + amp*np.sin(freq*2*np.pi*t)
    k2[lh_index] = lh0 + amp*np.sin(freq*2*np.pi*t)
    return tuple(k2)

# Create and initiate the integrator.
cmvi = PuppetMVI(puppet, collision_surfaces, puppet_strings, kin_func=h_func)
cmvi.initialize_state(0.0, puppet.q, np.zeros(cmvi.nd))

# Simulate the dynamics.
tf = 5.0
dt = 0.05
(t, q, lam) = tc.simulate(cmvi, tf, dt)

# Visualize the results.
pos = [12.8121004 , -10.37581321, 7.70056414] 
ang = [2.42159265, 0.06 , 0.]

cmvi.prepare_to_visualize()
if len(t)>1:
    trep.visual.visualize_3d([PuppetVisual(puppet, t, q)], camera_pos=pos, camera_ang=ang) 
else:
    print 'Displaying initial configuration only'
    trep.visual.visualize_3d([PuppetVisual(puppet, [0, 1], [puppet.q, puppet.q])], camera_pos=pos, camera_ang=ang) 
