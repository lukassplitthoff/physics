import numpy as np
import scipy.constants as pyc

from physics.mesoscopic.theorySC import Superconductivity


class JosephsonJunction(object):
    def __init__(self, delta, transmission, phasediff, eta):
        self.delta = delta
        self.transmission = transmission
        self.phasediff = phasediff
        self.eta = eta

    def calc_spectrum_short(self):
        return self.delta * np.sqrt(1 - self.transmission * (np.sin(self.phasediff / 2))**2)

    def calc_wavefun_short(self):
        return 1 / self.eta * np.sqrt(self.transmission) * np.abs(np.sin(self.phasediff / 2))

    def calc_supercurrent_short(self):
        return -2 * pyc.e / pyc.hbar * np.gradient(JosephsonJunction.calc_spectrum_short(self), self.phasediff)

    def calc_Josephson_inductance(self):
        return 1/np.gradient(JosephsonJunction.calc_supercurrent_short(self), self.phasediff)


def _example():

    sc = Superconductivity(delta=1, phi=np.pi/2)
    print('delta', sc.delta)
    print('phi', sc.phi)
    op = sc.calc_op_sc()
    print('order parameter', op)

    JJ = JosephsonJunction(delta=1, transmission=0.5, phasediff=np.linspace(0,2*np.pi, 1001), eta=1e-9)
    print('delta', JJ.delta)
    print('transmission', JJ.transmission)
    print('phasediff', JJ.phasediff)
    energy = JJ.calc_spectrum_short()
    print('energy', energy)
    wavefun = JJ.calc_wavefun_short()
    print('wavefunction', wavefun)
    supercurrent = JJ.calc_supercurrent_short()
    print('supercurrent', supercurrent)
    Josephson_inductance = JJ.calc_Josephson_inductance()
    print('Josephson_inductance', Josephson_inductance)


if __name__ == '__main__':
    _example()