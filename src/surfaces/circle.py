import trep
import numpy as np
from surface import Surface

# Import the visualization module. If this fails, the constraint will not
#   be drawn in the simulation.
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False


class CircleConstraint2(trep.Constraint):
    """ The actual constraint to be added once the mass is on the surface. """
    def __init__(self, system, Y, Z, R, frame):
        trep.Constraint.__init__(self, system, name='constraint_' + frame.name )
        self.Y = Y
        self.Z = Z
        self.R = R
        self.frame = frame

    def __repr__(self):
        return "<CircleConstraint %s Y=%f Z=%f R=%f>" %(
                self.frame.name, self.Y, self.Z, self.R)

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


class Circle2:
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


class Circle(Surface):
    def __init__(self, system, frame, Y=0, Z=0, R=1, name=None, tolerance=1e-10, ):
        Surface.__init__(self, system, name, tolerance)
        self.frame = frame
        self.Y = Y
        self.Z = Z
        self.R = R
    
    def phi(self):
        y, z = self.frame.p()[1:3]
        return (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def phi_dq(self, q_i):
        y, z = self.frame.p()[1:3]
        dy, dz = self.frame.p_dq(q_i)[1:3]
        return 2*(y-self.Y)*dy + 2*(z-self.Z)*dz 

    def phi_dqdq(self, q1, q2):
        y, z = self.frame.p()[1:3]
        dy1, dz1 = self.frame.p_dq(q1)[1:3]
        dy2, dz2 = self.frame.p_dq(q2)[1:3]
        dy12, dz12 = self.frame.p_dqdq(q1, q2)[1:3]
        return 2.0*(dy2*dy1+(y-self.Y)*dy12 + dz2*dz1 + (z-self.Z)*dz12)
        
    def phi_dqdqdq(self, q1, q2, q3):
        y, z = self.frame.p()[1:3]
        dy1, dz1 = self.frame.p_dq(q1)[1:3]
        dy2, dz2 = self.frame.p_dq(q2)[1:3]
        dy3, dz3 = self.frame.p_dq(q3)[1:3]
        dy23, dz23 = self.frame.p_dqdq(q2, q3)[1:3]
        dy13, dz13 = self.frame.p_dqdq(q1, q3)[1:3]
        dy12, dz12 = self.frame.p_dqdq(q1, q2)[1:3]
        dy123, dz123 = self.frame.p_dqdqdq(q1, q2, q3)[1:3]
        return 2.0*(dy23*dy1+dy1*dy13+dy3*dy12+(y-self.Y)*dy123+
                  dz23*dz1+dz2*dz13+dz3*dz12+(z-self.Z)*dz123) 

    def __repr__(self):
        return "<Circle Surface frame=%s>" %self.frame.name

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

