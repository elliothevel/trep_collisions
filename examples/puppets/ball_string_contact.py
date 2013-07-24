import trep
import trep_collisions as tc
from puppetmvi import *
from strings import *

s = trep.System()
s.import_frames([trep.tz(3.0, name='ref'),[trep.tz('z', mass=1.0, name='ball')]])
trep.potentials.Gravity(s)
c1 = trep.constraints.Distance(s, s.get_frame('ref'), s.get_frame('ball'), 1.0)
s1 = PuppetString(c1)

surface = tc.surfaces.Ground(s, dim=2)


pmvi = PuppetMVI(s, surface, puppet_strings=[s1], impact_frames=[s.get_frame('ball')])
pmvi.initialize_state(0.0, [1.0], [0.0])

tf = 1.0
dt = 0.1
(t, q, lam) = tc.simulate(pmvi, tf, dt)


print pmvi.system.constraints
pmvi.prepare_to_visualize()
trep.visual.visualize_3d([trep.visual.VisualItem3D(s, t, q)])
