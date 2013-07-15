import trep
import numpy as np

# Import the visualization module. If this fails, the constraint will not
#   be drawn in the simulation.
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False


class CircleConstraint(trep.Constraint):
    """ The actual constraint to be added once the mass is on the surface. """
    def __init__(self, system, Y, Z, R, frame):
        trep.Constraint.__init__(self, system, name='constraint_' + frame.name )
        self.Y = Y
        self.Z = Z
        self.R = R
        self.frame = frame

    def h(self):
        y, z = self.frame.p()[1:3]
        return (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def h_dq(self, q_i):
        y, z = self.frame.p()[1:3]
        dy, dz = self.frame.p_dq(q_i)[1:3]
        return 2*(y-self.Y)*dy + 2*(z-self.Z)*dz

    if _opengl:
        def opengl_draw(self):
            """ Draws a simple circle. """
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            glColor3f(1.0, 0.0, 0.0)
            glDisable(GL_LIGHTING)
            glLineWidth(2.0)
            glBegin(GL_LINE_LOOP)
            for i in range(360):
                angle = 2*np.pi*i/360
                glVertex3f(0.0,self.Y+self.R*np.cos(angle),self.Z+self.R*np.sin(angle))   
            glEnd()
            glPopAttrib()


class Circle:
    """ Defines a circular boundary in the plane. This function will give the value of
        phi without actually constraining the system. """    
    def __init__(self, system, Y, Z, R, tolerance=1e-10):
        self.Y = Y
        self.Z = Z
        self.R = R
        self.system = system
        self.tolerance = tolerance
        self.frame = None

    def set_frame(self, frame):
        """ Updates the frame that contacts the surface. """
        self.frame = self.system.get_frame(frame)	    

    def phi(self):
        y, z = self.frame.p()[1:3]
        return (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def phi_dq(self, q_i):
        y, z = self.frame.p()[1:3]
        dy, dz = self.frame.p_dq(q_i)[1:3]
        return 2*(y-self.Y)*dy + 2*(z-self.Z)*dz        
                
    def add_constraint_to_system(self, frame):
        """ Adds the constraint to the system. """
        CircleConstraint(self.system, self.Y, self.Z, self.R, frame)
