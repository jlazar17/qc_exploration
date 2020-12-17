from qiskit import Aer, execute, ClassicalRegister, QuantumRegister, QuantumCircuit
from collections.abc import Iterable
import numpy as np

from physicsconstants import OscParams

PMNS_param        = (-0.6031, 7.412,   0.7966,  1.0139, 0.7053, -8.065)
PMNS_dagger_param = (-0.7053, -1.3599, 0.7966, -1.0139, 0.6031, 2.0125)
op = OscParams()

def phi(m2, LoE):
    r'''
    params
    ______
    m2  (float): Mass squared difference [ev^2]
    LoE (float): Baseline over energy [km/GeV]
    
    returns
    _______
    phi (float): 
    '''

    coef = 2.534
    phi  = coef * m2 * LoE
    return phi

class ThreeNuOscillator:

    def __init__(self):
        self.qreg   = QuantumRegister(2)
        self.creg   = ClassicalRegister(2)
        self.qc     = QuantumCircuit(self.qreg, self.creg)
        self.counts = None
        

    def apply_PMNS(self, param=PMNS_param):
        r'''
    
        '''
        alpha, beta, gamma, delta, epsilon, zeta = param
        self.qc.u(epsilon, 0, 0, qreg[0])
        self.qc.u(zeta, 0, 0, qreg[1])
        self.qc.cnot(qreg[0], qreg[1])
        self.qc.u(gamma, 0, 0, qreg[0])
        self.qc.u(delta, 0, 0, qreg[1])
        self.qc.cnot(qreg[0], qreg[1])
        self.qc.u(alpha, 0, 0, qreg[0])
        self.qc.u(beta, 0, 0, qreg[1])

    def apply_PMNS_dagger(self, param=PMNS_dagger_param)
        r'''
    
        '''
        alpha, beta, gamma, delta, epsilon, zeta = param
        self.qc.u(epsilon, 0, 0, qreg[1])
        self.qc.u(zeta, 0, 0, qreg[0])
        self.qc.cnot(qreg[0], qreg[1])
        self.qc.u(gamma, 0, 0, qreg[0])
        self.qc.u(delta, 0, 0, qreg[1])
        self.qc.cnot(qreg[0], qreg[1])
        self.qc.u(alpha, 0, 0, qreg[0])
        self.qc.u(beta, 0, 0, qreg[1])


    def propoagate(self, LoE, m12=op.deltam12, m13=op.deltam13):
        r'''

        '''
        self.qc.rz(phi(m12, LoE), self.qreg[0])
        self.qc.rz(phi(m13, LoE), self.qreg[1])

if __name__=='__main__':
    loee = np.linsapce(0, 10000, 1000)
    
    n = 10000
    for i, LE in enumerate(loee):
    tno = ThreeNuOscillator()   

