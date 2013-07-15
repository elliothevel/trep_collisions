import trep

# Import the visualization module. If this fails, the constraint will not
#   be drawn in the simulation.
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
    _opengl = True
except:
    _opengl = False


class GroundConstraint(trep.constraints.PointOnPlane):
    """ Uses point on plane constraint but adds visualization. """
    def __init__(self, system, frame, dim, lims=(-5,5)):
        trep.constraints.PointOnPlane.__init__(self, system, system.world_frame, (0,0,-1), 
                                      frame, name='constraint_'+frame.name)
        self.lims = lims
        self.dim = dim

    if _opengl:
        def opengl_draw(self):
            """ Draws the line for the ground in 2-D or a plane for the ground in 3-D. """

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
    

class Ground:
    """ A common collision surface. This is a plane at z=0.  """ 
    def __init__(self, system, tolerance=1e-10, dim=2, lims=(-3,3)):
        self.system = system
        self.frame = None
        self.tolerance=tolerance
        self.dim = dim
        self.lims = lims

    def set_frame(self, frame):
        self.frame = self.system.get_frame(frame)
    
    def phi(self):
        return self.frame.p()[2]

    def phi_dq(self,q_i):
        return self.frame.p_dq(q_i)[2]

    def add_constraint_to_system(self, frame):
        """ Adds the constraint to the system. """
        GroundConstraint(self.system, frame, self.dim, lims=self.lims)
