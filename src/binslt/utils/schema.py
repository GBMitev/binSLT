from enum import Enum

class Schema(Enum):

    STATES = (
        "Columns for minimal ExoMol .states files\n"
        "NN = Counting Number\n"
        "E = Term Value (cm-1)\n"
        "gns = Nuclear Degeneracy\n"
        "J = Rotational Quantum Number\n"
        "tau = Parity +/-\n"
        "e/f = Rotationless Parity e/f\n"
        "Manifold = State label e.g. 12Pi\n"
        "Lambda = Molecular Orbital Angular Momentum\n"
        "Sigma = Molecular Spin Angular Momentum\n"
        "Omega = Total Angular Momentum\n"
        ,
        ["NN","E","gns","J","tau","e/f","Manifold","v","Lambda","Sigma","Omega"]
        )

    def describe(self):
        return self.value[0]
    
    def values(self):
        return self.value[1]

