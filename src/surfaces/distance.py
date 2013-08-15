import trep
from surface import Surface

try:
    from OpenGL.GLUT import *
    from OpenGL.GLU import *
    from OpenGL.GL import *
    _opengl = True
except:
    _opengl = False



class Distance(Surface):
    def __init__(self, system, frame1, frame2, dist, invalid='short'):
        Surface.__init__(self, system, 'Distance')
        self.frame1 = self.system.get_frame(frame1)
        self.frame2 = self.system.get_frame(frame2)
        self.tape_measure = trep.TapeMeasure(system, [frame1, frame2])

        # The invalid keyword controls whether phi < 0 when the distance
        # is less than or greater than the specified distance.
        assert invalid in ['short', 'long'], "Incorrect 'invalid' type"
        if invalid == 'short':
            self.sgn = 1.0
        else:
            self.sgn = -1.0

        if isinstance(dist, str):
            self.config = trep.Config(system=self.system, name=dist, kinematic=True)
        else:
            self.config = None
            self.dist = dist

        self.dist = dist
        self.sign = -self.sgn

    def __repr__(self):
        return "<Distance - frame1=%s, frame2=%s>" %(self.frame1.name, self.frame2.name)

    def measure_length(self):
        return self.tape_measure.length()

    def phi(self):
        if self.config:
            return self.sgn*(self.tape_measure.length() - self.config.q)  
        else:    
            return self.sgn*(self.tape_measure.length() - self.dist)

    def phi_dq(self, q1):  
        if (self.config is None) or self.config != q1:
            return self.sgn*self.tape_measure.length_dq(q1)
        return self.sgn*(self.tape_measure.length_dq(q1)-1)

    def phi_dqdq(self, q1, q2):
        return self.sgn*self.tape_measure.length_dqdq(q1, q2)

    if _opengl:
        def opengl_draw(self):
            frame1 = self.frame1.g()
            frame2 = self.frame2.g()

            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT )
            glColor3f(0.0, 0.0, 1.0)
            glDisable(GL_LIGHTING)
            glBegin(GL_LINES)
            glVertex3f(frame1[0][3], frame1[1][3], frame1[2][3])
            glVertex3f(frame2[0][3], frame2[1][3], frame2[2][3])    
            glEnd()
            glPopAttrib() 
