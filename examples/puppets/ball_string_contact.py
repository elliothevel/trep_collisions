import trep
import trep_collisions as tc
from strings2 import *

s = trep.System()
s.import_frames([trep.tz('z_ref', name='ref', kinematic=True),[trep.tz('z', mass=1.0, name='ball')]])
trep.potentials.Gravity(s)

surfaces = tc.surfaces.global_surfaces(tc.surfaces.Ground, s, kwargs={'dim': 2}, impact_frames=[s.get_frame('ball')])
#surfaces += [PuppetString(s, s.get_frame('ref'), s.get_frame('ball'), 1.0)]

trep.potentials.LinearSpring(s, 'ball', 'ref', k=1000.0, x0=1.0)
trep.forces.Damping(s, 0.1, {'z': 100.0})

def kin_func(t):
    if t <= 2.5:
        z = 3.0-1.0*t
    else:
        z = -2.0+1.0*t
    return (z, )


pmvi = tc.CollisionMVI(s, surfaces, kin_func=kin_func)
s.get_config('z').q = - 1.0
s.get_config('z_ref').q = 3.0
pmvi.initialize_state(0.0, s.q, [0.0])

tf = 6.0
dt = 0.1
(t, q, lam) = tc.simulate(pmvi, tf, dt)

pmvi.prepare_to_visualize()
trep.visual.visualize_3d([trep.visual.VisualItem3D(s, t, q)])
