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

class FixedTrajectory(trep.Constraint):
    def __init__(self, system, name, amp, freq):
        trep.Constraint.__init__(self, system)
        self.config = system.get_config(name)
        self.amp = amp
        self.freq = freq
        self.offset = self.config.q

    def h(self):
        return (self.config.q - self.amp*np.sin(2*np.pi*self.freq*self.system.t)-self.offset)
    
    def h_dq(self,q_i):
        if self.config == q_i:
            return 1.0
        return 0.0

    def h_dqdq(self,q_i,q_j):
        return 0.0

    #def h_dt(self):
    #    return -self.amp*2*np.pi*self.freq*np.cos(2*np.pi*self.freq*self.system.t)

amp, freq = -2., 0.2
FixedTrajectory(puppet, 'Frame Z', amp, freq) 
puppet.satisfy_constraints()



# Initiate the integrator.
cmvi = PuppetMVI(puppet, surface, impact_frames=puppet.feet_frames[:2], puppet_strings=puppet.strings)
#puppet.get_config('RKneeTheta').q += 0.01
cmvi.initialize_state(0.0, puppet.q, np.zeros(cmvi.nd))

# Define a function that will move the puppet up and down. This is keeping the
# strings a constant length and moving the overhead frame up and down.



amp, freq = -2., 0.1
h0 = puppet.get_config('Frame Z').q
def h_func(q, t):
    h = h0 + amp*np.sin(freq*2*np.pi*t)
    return (h, )

# Simulate the dynamics.
tf = 7.0
dt = 0.01

print '\n\n'
(t, q, lam) = tc.simulate(cmvi, tf, dt)#, h_func)

# Visualize the results.
pos = [12.8121004 , -10.37581321, 7.70056414] 
ang = [2.42159265, 0.06 , 0.]
cmvi.prepare_to_visualize()
trep.visual.visualize_3d([PuppetVisual(puppet, t, q)], camera_pos=pos, camera_ang=ang)    
