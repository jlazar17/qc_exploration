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

    def __init__(self, init_state):
        self.qreg   = QuantumRegister(2)
        self.creg   = ClassicalRegister(2)
        self.qc     = QuantumCircuit(self.qreg, self.creg)
        if init_state=='nue':
            pass
        elif init_state=='numu':
            self.qc.x(self.qreg[1])
        elif init_state=='nutau':
            self.qc.x(self.qreg[0])
        else:
            print('init_state %s not recognized. Please reinitialize' % init_state)
        self.counts = None
        

    def apply_PMNS(self, param=PMNS_param):
        r'''
    
        '''
        alpha, beta, gamma, delta, epsilon, zeta = param
        self.qc.u(beta, 0, 0, self.qreg[1])
        self.qc.u(alpha, 0, 0, self.qreg[0])
        self.qc.cnot(self.qreg[1], self.qreg[0])
        self.qc.u(delta, 0, 0, self.qreg[1])
        self.qc.u(gamma, 0, 0, self.qreg[0])
        self.qc.cnot(self.qreg[1], self.qreg[0])
        self.qc.u(zeta, 0, 0, self.qreg[1])
        self.qc.u(epsilon, 0, 0, self.qreg[0])

    def apply_PMNS_dagger(self, param=PMNS_dagger_param):
        r'''
    
        '''
        alpha, beta, gamma, delta, epsilon, zeta = param
        self.qc.u(beta, 0, 0, self.qreg[1])
        self.qc.u(alpha, 0, 0, self.qreg[0])
        self.qc.cnot(self.qreg[1], self.qreg[0])
        self.qc.u(delta, 0, 0, self.qreg[1])
        self.qc.u(gamma, 0, 0, self.qreg[0])
        self.qc.cnot(self.qreg[1], self.qreg[0])
        self.qc.u(zeta, 0, 0, self.qreg[1])
        self.qc.u(epsilon, 0, 0, self.qreg[0])


    def propoagate(self, LoE, m12=op.deltam12, m13=op.deltam3l):
        r'''

        '''
        self.qc.rz(phi(m12, LoE), self.qreg[1])
        self.qc.rz(phi(m13, LoE), self.qreg[0])

    def measure(self):
        self.qc.measure(self.qreg, self.creg)

if __name__=='__main__':
    loee = np.linspace(0, 1200, 21)
    
    n = 10000
    results = np.zeros((4,len(loee)))
    for i, LE in enumerate(loee):
        tno = ThreeNuOscillator('numu')
        tno.apply_PMNS_dagger()
        tno.propoagate(LE)
        tno.apply_PMNS()
        tno.measure()
        job = execute(tno.qc, Aer.get_backend('qasm_simulator'), shots=n)
        counts = job.result().get_counts(tno.qc)
        for j, (key, val) in enumerate(sorted(counts.items())):
            results[j, i] = float(val)/n
    np.save('three_neutrino', results)
