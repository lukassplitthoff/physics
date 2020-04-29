import numpy as np

from mesoscopic.theorySC import Superconductivity


class JosephsonJunction(object):
    def __init__(self, delta, transmission, phasediff):
        self.delta = delta
        self.transmission = transmission
        self.phasediff = phasediff

    def calc_spectrum_short(self):
        return self.delta*np.sqrt(1-self.transmission*(np.sin(self.phasediff/2))**2)


def _example():

    sc = Superconductivity(delta=1, phi=np.pi/2)
    print('delta', sc.delta)
    print('phi', sc.phi)
    op = sc.calc_op_sc()
    print('order parameter', op)

    JJ = JosephsonJunction(delta=1, transmission=0.5, phasediff=np.pi)
    print('delta', JJ.delta)
    print('transmission', JJ.transmission)
    print('phasediff', JJ.phasediff)
    energy = JJ.calc_spectrum_short()
    print('energy', energy)


if __name__ == '__main__':
    _example()