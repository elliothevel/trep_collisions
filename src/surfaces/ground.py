import trep
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
  

class Ground(Surface):
    def __init__(self, system, frame, name=None, tolerance=1e-10, dim=2, lims=(-3,3)):
        Surface.__init__(self, system, name, tolerance)
        self.frame = frame
        self.dim = dim
        self.lims = lims
    
    def phi(self):
        return self.frame.p()[2]

    def phi_dq(self, q1):
        return self.frame.p_dq(q1)[2]

    def phi_dqdq(self, q1, q2):
        return self.frame.p_dqdq(q1, q2)[2]

    def phi_dqdqdq(self, q1, q2, q3):
        return self.frame.p_dqdqdq(q1, q2, q3)

    def __repr__(self):
        return "<Ground Surface frame=%s>" %self.frame.name

    if _opengl:
        def opengl_draw(self):
            if self.dim == 2:
                glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
                glColor3f(1.0, 0.0, 0.0)
                glDisable(GL_LIGHTING)
                glLineWidth(2.0)

                glBegin(GL_LINES)
                glVertex3f(0.0, self.lims[0], 0.0)
                glVertex3f(0.0, self.lims[1], 0.0)    
                glEnd()

                glPopAttrib()

            if self.dim == 3:  
                x_min, x_max = self.lims
                y_min, y_max = self.lims
                glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
                glColor3f(96/255.0, 96/255.0, 96/255.0)
                glDisable(GL_LIGHTING)
                glLineWidth(2.0)

                glBegin(GL_POLYGON)
                glVertex3f(x_min, y_min, 0.0)
                glVertex3f(x_min, y_max, 0.0)
                glVertex3f(x_max, y_max, 0.0)
                glVertex3f(x_max, y_min, 0.0)
                glEnd()

                glPopAttrib() 
