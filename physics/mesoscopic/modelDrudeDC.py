import numpy as np
import scipy.constants as pyc


class DrudeModelDC(object):
    def __init__(self, density, meanfreetime, charge=pyc.e, mass=pyc.m_e, frequency=0.0):
        self.density = density
        self.charge = charge
        self.tau = meanfreetime
        self.mass = mass
        self.frequency_angular = 2 * np.pi * frequency

    def calc_complexconductivity(self):
        return self.density * self.charge ** 2 * self.tau / (self.mass * (1 + 1j * self.frequency_angular * self.tau))

    def calc_currentdensity(self, electricfield):
        return DrudeModelDC.calc_complexconductivity(self) * electricfield

    def calc_kineticinductance(self, length, crosssection):
        return self.mass / (self.density * self.charge**2) * length / crosssection


def _example():

    wire = DrudeModelDC(density=1e25, meanfreetime=1e-14, charge=pyc.e, mass=pyc.m_e, frequency=4.0e9)
    print('complex conductivity', wire.calc_complexconductivity())
    print('current density', wire.calc_currentdensity(electricfield=2))
    print('kinetic inductance', wire.calc_kineticinductance(length=2e-6, crosssection=np.pi * 40e-9 ** 2))


if __name__ == '__main__':
    _example()



