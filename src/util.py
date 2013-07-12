import numpy as np


def set_system_state(mvi, i):
    """ Sets the state of sytem to either 1 or 2. """
    if i==1:
        mvi.system.set_state(q=mvi.q1, t=mvi.t1)
    if i==2:
        mvi.system.set_state(q=mvi.q2, t=mvi.t2)


def D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([0.5*mvi.system.L_dq(q)*dt-mvi.system.L_ddq(q)
                     for q in mvi.system.dyn_configs])


def D2L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([0.5*mvi.system.L_dq(q)*dt+mvi.system.L_ddq(q)
                    for q in mvi.system.dyn_configs])


def D2D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[0.25*dt*mvi.system.L_dqdq(q1,q2)
                     +0.5*mvi.system.L_ddqdq(q1,q2)
                     -0.5*mvi.system.L_ddqdq(q2,q1)
                     -1.0/dt*mvi.system.L_ddqddq(q1,q2)
                    for q1 in mvi.system.dyn_configs]
                    for q2 in mvi.system.dyn_configs])


def fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([sum([force.f(q)*dt
                   for force in mvi.system.forces])
                   for q in mvi.system.dyn_configs])


def D2fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[sum([0.5*force.f_dq(q1,q2)*dt
                    +force.f_ddq(q1,q2)
                    for force in mvi.system.forces])
                    for q2 in mvi.system.dyn_configs]
                    for q1 in mvi.system.dyn_configs])


def D4D1L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[0.5*mvi.system.L_dq(q1)
                      -0.5*np.dot(mvi.system.dq, [mvi.system.L_ddqdq(q2,q1)
                                 for q2 in mvi.system.configs])
                      +1.0/dt*np.dot(mvi.system.dq,[mvi.system.L_ddqddq(q1,q2)
                                 for q2 in mvi.system.configs])]
                     for q1 in mvi.system.dyn_configs])

    
def D4fm2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([[sum([force.f(q1) for force in mvi.system.forces])
                      -np.dot(mvi.system.dq, [sum([force.f_dq(q1,q2)
                    for force in mvi.system.forces])
                    for q2 in mvi.system.configs])]
                    for q1 in mvi.system.dyn_configs])


def D4L2(mvi):
    return mvi.system.L()-np.dot(mvi.system.dq,
            [mvi.system.L_ddq(q) for q in mvi.system.configs])


def D3L2(mvi):
    return -mvi.system.L()+np.dot(mvi.system.dq,
            [mvi.system.L_ddq(q) for q in mvi.system.configs])


def D2D3L2(mvi):
    dt = mvi.t2-mvi.t1
    return np.array([-0.5*mvi.system.L_dq(q1)
                    +np.dot(mvi.system.dq, [0.5*mvi.system.L_ddqdq(q1,q2)
		            +1.0/dt*mvi.system.L_ddqddq(q1,q2) for q2 in mvi.system.configs])
                    for q1 in mvi.system.dyn_configs])


def h(mvi):
    return np.array([c.h() for c in mvi.system.constraints])


def Dh(mvi):
    # This needs a transpose to handle the case of no constraints,
    # where the shapes of the empty arrays need to match.
    return np.transpose(np.array([[c.h_dq(q) for c in mvi.system.constraints]
                                   for q in mvi.system.dyn_configs]))

    
def phi(mvi):
    mvi.surface.system.set_state(q=mvi.system.q, t=mvi.system.t)
    return mvi.surface.phi()


def Dphi(mvi, q_i=None):
    mvi.surface.system.set_state(q=mvi.system.q, t=mvi.system.t)
    if q_i is not None:
        return mvi.surface.h_dq(q_i)
    return np.array([[mvi.surface.phi_dq(q)] for q in mvi.system.dyn_configs])


def D2h(mvi):
    mvi.surface.system.set_state(q=mvi.system.q, t=mvi.system.t)
    Dh = np.zeros(mvi.nc)
    for i in range(mvi.nc):
        constraint = mvi.system.get_constraint(i)
        try:
            Dh[i] = constraint.h_dt()
        except:
            pass
    return Dh
    
    
def save_state(mvi):
    """ Saves the integrator's states in a dictionary. """
    state = {'q0': mvi.q0, 'q1': mvi.q1, 'q2': mvi.q2,
             'p0': mvi.p0, 'p1': mvi.p1, 'p2': mvi.p2,
             't0': mvi.t0, 't1': mvi.t1, 't2': mvi.t2,
             'lambda0': mvi.lambda0, 'lambda1': mvi.lambda1}
    return state


def restore_state(mvi, state):
    """ Resets the integrator's state from a dictionary. """
    mvi.t1 = state['t1']; mvi.t2 = state['t2']
    mvi.q1 = state['q1']; mvi.q2 = state['q2']
    mvi.p1 = state['p1']; mvi.p2 = state['p2']        