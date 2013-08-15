import trep
import trep_collisions as tc
from puppetmvi import PuppetMVI
import numpy as np

# Creat the system and add the frames.
s = trep.System()
s.import_frames([trep.tz('z_ref', name='ref', kinematic=True),[trep.tz('z', mass=1.0, name='ball')]])
trep.potentials.Gravity(s)

# Create the surfaces - the ground and a string connecting the two frames.
collision_surfaces = tc.surfaces.global_surfaces(tc.surfaces.Ground, s, kwargs={'dim': 2}, impact_frames=[s.get_frame('ball')])
strings = [tc.surfaces.Distance(s, 'ref', 'ball', 1.0, invalid='short')]

# Function for kinematically controlling the hand.
def kin_func(t):
    if t <= 2.5:
        z = 3.0-1.0*t
    else:
        z = -2.0+1.0*t
    return (z, )

pmvi = PuppetMVI(s, collision_surfaces, strings, kin_func=kin_func)
s.get_config('z').q = - 1.0
s.get_config('z_ref').q = 3.0
pmvi.initialize_state(0.0, s.q, [0.0], lambda1=np.zeros(pmvi.nc))
pmvi.add_all_strings()

tf = 6.0
dt = 0.013
(t, q, lam) = tc.simulate(pmvi, tf, dt)

pmvi.prepare_to_visualize()
print s.constraints
trep.visual.visualize_3d([trep.visual.VisualItem3D(s, t, q)])
