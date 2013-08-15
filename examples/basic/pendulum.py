"""
This example solves the problem of a multilink pendulum impacting a circle in 2D.
It's goal is to demonstrate how to use the trep_collisions package.
"""

import trep
from trep.visual import *
import trep_collisions as tc

# Define the system in trep. There is gravity and damping at each of the links.
L = .50
M = 1.0
s = trep.System()
frames = ([trep.rx('theta1'),[
           trep.tz(-L, mass=M), [trep.rx('theta2'), [
           trep.tz(-L, mass=M), [trep.rx('theta3'), [
           trep.tz(-L, mass=M), [trep.rx('theta4'), [
           trep.tz(-L, mass=M), [trep.rx('theta5'), [
           trep.tz(-L, mass=M), [trep.rx('theta6'), [
           trep.tz(-.50, mass=1.0)]]]]]]]]]]]])
s.import_frames(frames)
trep.potentials.Gravity(s)
trep.forces.Damping(s, .50)

# Define the collision surface.
surface = tc.surfaces.global_surfaces(tc.surfaces.Circle, s, kwargs={'Y': -2, 'Z': -1.7, 'R': 1.0}, impact_frames=s.masses)

# Create the variational integrator with collisions.                   
mvi = tc.CollisionMVI(s, surface, release_method='interp', impact_method='root')
#mvi._releases_off = True

# Give the system some in initial conditions.
s.get_config('theta1').q = -2.0
mvi.initialize_state(0.0, s.q, [0.0]*s.nQd)

# Simulate the dynamics.
tf = 1.
dt = 0.02
t, q, lam= tc.simulate(mvi, tf, dt)

# Visualize.
mvi.prepare_to_visualize()
visualize_3d([VisualItem3D(s, t, q)])
