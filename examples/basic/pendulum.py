"""
This example solves the problem of a multilink pendulum impacting a circle in 2D.
It's goal is to demonstrate how to use the trep_collisions package.
"""

import trep
from trep.visual import *
import trep_collisions as tc

# Define the system in trep. There is gravity and damping at each of the links.
s = trep.System()
frames = ([trep.rx('theta1'),[
           trep.tz(-.50, mass=1.0), [trep.rx('theta2'), [
           trep.tz(-.50, mass=1.0), [trep.rx('theta3'), [
           trep.tz(-.50, mass=1.0), [trep.rx('theta4'), [
           trep.tz(-.50, mass=1.0), [trep.rx('theta5'), [
           trep.tz(-.50, mass=1.0), [trep.rx('theta6'), [
           trep.tz(-.50, mass=1.0)]]]]]]]]]]]])
s.import_frames(frames)
trep.potentials.Gravity(s)
trep.forces.Damping(s, 1.0)

# Define the collision surface.
surface = tc.surfaces.global_surfaces(tc.surfaces.Circle, s, kwargs={'Y': -2, 'Z': -1.7, 'R': 1})

# Create the variational integrator with collisions.                   
mvi = tc.CollisionMVI(s, surface, release_method='interp')

# Give the system some in initial conditions.
s.get_config('theta1').q = 2.0
mvi.initialize_state(0.0, s.q, [0.0]*s.nQd)

# Simulate the dynamics.
tf = 6.0
dt = 0.01
t, q, lam= tc.simulate(mvi, tf, dt)

# Visualize.
mvi.prepare_to_visualize()
visualize_3d([VisualItem3D(s, t,q)])
