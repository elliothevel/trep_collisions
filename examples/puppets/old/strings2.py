import trep
import trep_collisions as tc


try:
    from OpenGL.GLUT import *
    from OpenGL.GLU import *
    from OpenGL.GL import *
    _opengl = True
except:
    _opengl = False


class PuppetString2(trep.constraints.Distance):
    def __init__(self, system, frame1, frame2, distance, tolerance=1e-10):
        system.hold_structure_changes()
        super(PuppetString, self).__init__(system, frame1, frame2, distance) 
        self.system._constraints = self.system._constraints[:-1]
        self.system.resume_structure_changes()
        self.active = False
        
    def deactivate_constraint(self):
        self.system._constraints = tuple([con for con in self.system.constraints if (con != self)]) # TUPLE!!
        self.active = False 

    def activate_constraint(self):
        self.system._add_constraint(self)
        self.active = True

    def phi(self):
        return self.h()

    def phi_dq(self, q1):
        return self.h_dq(q1)

    #def h(self):
    #    return self.phi()

    #def h_dq(self, q1):
    #    return self.phi_dq(q1)
    

class PuppetString(tc.surfaces.Surface):
    def __init__(self, system, frame1, frame2, dist):
        tc.surfaces.Surface.__init__(self, system, 'PuppetString')
        self.frame1 = self.system.get_frame(frame1)
        self.frame2 = self.system.get_frame(frame2)
        self.tape_measure = trep.TapeMeasure(system, [frame1, frame2])

        if isinstance(dist, str):
            self.config = trep.Config(system=self.system, name=dist, kinematic=True)
        else:
            self.config = None
            self.dist = dist

        self.dist = dist

    def __repr__(self):
        return "<PuppetString - frame1=%s, frame2=%s>" %(self.frame1.name, self.frame2.name)

    # We want phi < 0 when the length is greater than the given length.
    def phi(self):
        if self.config:
            return self.config.q - self.tape_measure.length()   
        else:    
            return self.dist - self.tape_measure.length()

    def phi_dq(self, q1):    
        return -self.tape_measure.length_dq(q1)

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

