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

class CosineConstraint(trep.Constraint):
    """  The actual constraint to be added once the mass is on the surface.  """
    def __init__(self, system, amp, freq, frame, ymin=-4*np.pi,ymax=4*np.pi):
        trep.Constraint.__init__(self, system, name='constraint_' + frame.name)
        self.ylims = [ymin, ymax]
        self.amp = amp
        self.freq = freq
        self.frame = frame

    def __repr__(self):
        return "<CosineConstraint %s amp=%f freq=%f>" %(
                self.frame.name, self.amp, self.freq)

    def h(self):
        y, z = self.frame.p()[1:3]
        return z-self.amp*np.cos(self.freq*y)

    def h_dq(self, q_i):
        (dydq,dzdq) = self.frame.p_dq(q_i)[1:3]  
        dphidy = self.amp*self.freq*np.sin(self.freq*self.frame.p()[1])   
        return dzdq + dphidy*dydq

    if _opengl:
        def opengl_draw(self):
            Y = np.linspace(self.ylims[0],self.ylims[1],100)
            Z = self.amp*np.cos(self.freq*Y)
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_LINE_BIT)
            glColor3f(1.0, 0.0, 0.0)
            glDisable(GL_LIGHTING)
            glLineWidth(2.0)
            glBegin(GL_LINE_STRIP)
            for i in range(len(Y)):
                glVertex3f(0.0,Y[i],Z[i])   
            glEnd()
            glPopAttrib()
   

class Cosine:
    """ Defines a cosine curve as a collision surface. phi>0 when z>amp*cos(freq*y). """ 
    def __init__(self, system, amp, freq, tolerance=1e-10):
        self.system = system
        self.amp = amp
        self.freq = freq
        self.frame = None
        self.tolerance = tolerance

    def __repr__(self):
        return '<CosineSurface freq=%s, amp=%s>' %(self.freq, self.amp)    

    def set_frame(self,node_name):
        self.frame = self.system.get_frame(node_name)
    
    def phi(self):
        y, z = self.frame.p()[1:3]
        return z-self.amp*np.cos(self.freq*y)

    def phi_dq(self, q_i):
        dydq, dzdq = self.frame.p_dq(q_i)[1:3]  
        dphidy = self.amp*self.freq*np.sin(self.freq*self.frame.p()[1])   
        return dzdq + dphidy*dydq

    def add_constraint_to_system(self, frame):
        """ Adds the constraint to the system. """
        CosineConstraint(self.system, self.amp, self.freq, frame)
