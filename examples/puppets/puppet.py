import trep
from trep import tx, ty, tz, rx, ry, rz
import trep.visual as visual
import trep_collisions as tc
import numpy as np


PUPPET_FRAMES = [
    tx('TorsoX'), [ty('TorsoY'), [tz('TorsoZ'), [
        rz('TorsoPsi'), [ry('TorsoTheta'), [rx('TorsoPhi',name='Torso'), [
            tz(-1.5, mass=50),
            tx(-1.011), [tz(0.658, name='Right Torso Hook')],
            tx( 1.011), [tz(0.658, name= 'Left Torso Hook')],
            tz(0.9, name='Head'), [tz(0.5, mass=(10,1,1,1))],
            # Define the left arm
            tx(1.3), [tz(0.4), [
                rz('LShoulderPsi'), [ry('LShoulderTheta'), [rx('LShoulderPhi', name='Left Shoulder'), [
                    tz(-0.95, name='Left Humerus', mass=(5,1,1,1)),
                    tz(-1.9), [
                        rx('LElbowTheta', name='Left Elbow'), [
                            tz(-1, name='Left Radius', mass=(4,1,1,1)),
                            tz(-2.001), [tx(0.14), [ty(-0.173, name='Left Finger')]]]]]]]]],
            # Define the right arm
            tx(-1.3), [tz(0.4), [
                rz('RShoulderPsi'), [ry('RShoulderTheta'), [rx('RShoulderPhi', name='Right Shoulder'), [
                    tz(-0.95, name='Right Humerus', mass=(5,1,1,1)),
                    tz(-1.9), [
                        rx('RElbowTheta', name='Right Elbow'), [
                            tz(-1, name='Right Radius', mass=(4,1,1,1)),
                            tz(-2.001), [tx(-0.14), [ty(-0.173, name='Right Finger')]]]]]]]]],
            # Define the left leg
            tx(0.5), [tz(-3.0), [
                rz('LHipPsi'), [ry('LHipTheta'), [rx('LHipPhi', name='Left Hip'), [
                    tz(-1.5, name='Left Femur', mass=(5,1,1,1)),
                    tz(-2.59), [ty(-0.322, name='Left Knee Hook')],
                    tz(-3.0), [
                        rx('LKneeTheta', name='Left Knee'), [
                            tz(-1.5, name='Left Tibia', mass=(4,1,1,1)),[tz(-1.5, name='Left Foot', mass=1.0),[ty(-1.0, name='Left Toe', mass=.50)]]]]]]]]],
            # Define the right leg
            tx(-0.5), [tz(-3.0), [
                rz('RHipPsi'), [ry('RHipTheta'), [rx('RHipPhi', name='Right Hip'), [
                    tz(-1.5, name='Right Femur', mass=(5,1,1,1)),
                    tz(-2.59), [ty(-0.322, name='Right Knee Hook')],
                    tz(-3.0), [
                        rx('RKneeTheta', name='right Knee'), [
                            tz(-1.5, name='Right Tibia', mass=(4,1,1,1)),[tz(-1.5, name='Right Foot', mass=1.0),[ty(-1.0, name='Right Toe', mass=0.50)]]]]]]]]],
          ]]]]]],  # End of puppet definition
    
    # Define the coordinate frames for the truck that the puppet is suspended from
    tz(24.0, name='Frame Plane'), [
        tx( 1, name='Left Torso Spindle'),
        tx(-1, name='Right Torso Spindle'),
        tx( 1.44), [ty(-2, name='Left Arm Spindle')],
        tx(-1.44), [ty(-2, name='Right Arm Spindle')],
        tx( 0.5),  [ty(-2.6, name='Left Leg Spindle')],
        tx(-0.5),  [ty(-2.6, name='Right Leg Spindle')]]
    ]


class Puppet(trep.System):
    def __init__(self, add_strings):
        trep.System.__init__(self)

        self.import_frames(PUPPET_FRAMES)
        trep.potentials.Gravity(self)
        trep.forces.Damping(self, 0.01)

        self.add_string_constraints(**add_strings)

        self.feet_frames = [self.get_frame('Right Foot'), self.get_frame('Left Foot'),
                            self.get_frame('Right Toe'),  self.get_frame('Left Toe')]

        self.init_config(add_strings)    
              

    def add_string_constraints(self, shoulders=True, hands=False, knees=False):
        """ Adds string constraints to the puppet. """

        strings = []
                
        if shoulders:  
            strings.extend([tc.surfaces.Distance(self, 'Left Torso Hook', 'Left Torso Spindle', 'LeftShoulderString', invalid='long'),
                            tc.surfaces.Distance(self, 'Right Torso Hook', 'Right Torso Spindle', 'RightShoulderString', invalid='long')])

        if hands: 
            strings.extend([tc.surfaces.Distance(self, 'Left Finger', 'Left Arm Spindle', 'LeftHandString', invalid='long'),
                            tc.surfaces.Distance(self, 'Right Finger', 'Right Arm Spindle', 'RightHandString', invalid='long')])


        if knees:   
            strings.extend([tc.surfaces.Distance(self, 'Left Knee Hook', 'Left Leg Spindle', 'LeftKneeString', invalid='long'),
                            tc.surfaces.Distance(self, 'Right Knee Hook', 'Right Leg Spindle', 'RightKneeString', invalid='long')])

        self.strings = strings
        self._structure_changed()


    def init_config(self, add_strings):
        """ Initialzes the configuration based on number of strings. Activates the string constraints. """

        q0 = {'TorsoZ': 10.0}
                
        if add_strings['shoulders']:        
            q0.update({'LeftShoulderString': 13.34, 'RightShoulderString': 13.34})

        # If there are strings on the hands, we want the elbows bent. Otherwise, straight at the sides.
        if add_strings['hands']:
            q0.update({'LElbowTheta' : -np.pi/2, 'RElbowTheta' : -np.pi/2,
                       'LeftHandString': 15.33, 'RightHandString': 15.33})
            
        self.q = q0

        self.hold_structure_changes()
        for string in self.strings:
            string.activate_constraint()
        self.resume_structure_changes()    
        self.satisfy_constraints(keep_kinematic=True)                    


# Class to visualize the puppet.
class PuppetVisual(visual.VisualItem3D):
    def __init__(self, *args, **kwds):
        super(PuppetVisual, self).__init__(*args, **kwds)

        self.attachDrawing('Torso',          visual.stlmodel('./puppet-stl/torso.stl').draw)
        self.attachDrawing('Left Shoulder',  visual.stlmodel('./puppet-stl/lefthumerus.stl').draw)
        self.attachDrawing('Right Shoulder', visual.stlmodel('./puppet-stl/righthumerus.stl').draw)
        self.attachDrawing('Left Elbow',     visual.stlmodel('./puppet-stl/leftradius.stl').draw)
        self.attachDrawing('Right Elbow',    visual.stlmodel('./puppet-stl/rightradius.stl').draw)
        self.attachDrawing('Right Hip',      visual.stlmodel('./puppet-stl/femur.stl').draw)
        self.attachDrawing('Left Hip',       visual.stlmodel('./puppet-stl/femur.stl').draw)
        self.attachDrawing('right Knee',     visual.stlmodel('./puppet-stl/tibia.stl').draw)
        self.attachDrawing('Left Knee',      visual.stlmodel('./puppet-stl/tibia.stl').draw)
        self.attachDrawing('Head',           visual.stlmodel('./puppet-stl/head.stl').draw)
        self.attachDrawing(None, self.draw_constraints)

    def draw_constraints(self):
        for x in self.system.constraints:
            x.opengl_draw()
