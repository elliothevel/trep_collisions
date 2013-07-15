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


class SphereConstraint(trep.Constraint):
    """ The actual constraint to be added once the mass is on the surface. """
    def __init__(self, system, X, Y, Z, R, frame):
        trep.Constraint.__init__(self, system, name='constraint_' + frame.name)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.R = R
        self.frame = frame

    def h(self):
        x, y, z = self.frame.p()
        return (x-self.X)**2 + (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def h_dq(self,q_i):
        x, y, z = self.frame.p() 
        dx, dy, dz = self.frame.p_dq(q_i)
        return 2*(x-self.X)*dx + 2*(y-self.Y)*dy + 2*(z-self.Z)*dz

    if _opengl:
        def opengl_draw(self):
            """ Draw a simple sphere. """
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            subdivisions = 50
            color = [1,0,0,1]
            quad = gluNewQuadric()
            gluQuadricNormals(quad, GLU_SMOOTH)
            glPushMatrix()
            glTranslatef(self.X, self.Y, self.Z)
            glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color)
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, .1)
            gluSphere(quad, self.R, subdivisions, subdivisions)
            glPopMatrix()
            glPopAttrib()
            gluDeleteQuadric(quad)


class Sphere:
    """ Defines a spherical boundary in 3D. This function will give the value of
        phi without actually constraining the system. """    
    def __init__(self, X, Y, Z, R, system, tolerance=1e-10):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.R = R
        self.system = system
        self.tolerance = tolerance
        self.frame = None 

    def set_frame(self, frame):
        """ Updates the node and frame that contacts the surface. """
        self.frame = self.system.get_frame(frame)

    def phi(self):
        x, y, z = self.frame.p()
        return (x-self.X)**2 + (y-self.Y)**2 + (z-self.Z)**2 - self.R**2

    def phi_dq(self,q_i):
        x, y, z = self.frame.p()
        dx, dy, dz = self.frame.p_dq(q_i)
        return 2*(x-self.X)*dx + 2*(y-self.Y)*dy + 2*(z-self.Z)*dz    
                
    def add_constraint_to_system(self, frame):
        """ Adds the constraint to the system. """
        SphereConstraint(self.system, self.X, self.Y, self.Z, self.R, frame)
    
