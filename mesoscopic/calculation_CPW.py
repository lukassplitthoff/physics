import numpy as np
import scipy.constants as pyc


class CPW(object):
    """capacitance per unit length in zero order approximation
    http://qucs.sourceforge.net/tech/node86.html

    Inductances and attenuation constant for a thin-film superconducting coplanar waveguide resonator
    http://dx.doi.org/10.1063/1.4773070"""
    def __init__(self, width, spacing, epsilon_r, L_kin_sq):
        self.width = width
        self.spacing = spacing
        self.epsilon_r = epsilon_r
        self.L_kin_sq = L_kin_sq

    def k(self):
        return self.width / (self.width + 2 * self.spacing)

    def calc_capacitance_air(self):
        if CPW.k(self) <= 1 / np.sqrt(2):
            kprime = np.sqrt(1 - CPW.k(self) ** 2)
            F = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
        else:
            F = np.log(2 * (1 + np.sqrt(CPW.k(self))) / (1 - np.sqrt(CPW.k(self)))) / np.pi

        return 2 * pyc.epsilon_0 * F

    def calc_capacitance_dielectric(self):

        if CPW.k(self) <= 1 / np.sqrt(2):
            kprime = np.sqrt(1 - CPW.k(self) ** 2)
            F = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
        else:
            F = np.log(2 * (1 + np.sqrt(CPW.k(self))) / (1 - np.sqrt(CPW.k(self)))) / np.pi

        return 2 * pyc.epsilon_0 * self.epsilon_r * F

    def calc_capacitance_total(self):
        return CPW.calc_capacitance_dielectric(self) + CPW.calc_capacitance_air(self)

    def calc_impedance_characteristic(self):

        if CPW.k(self) <= 1 / np.sqrt(2):
            kprime = np.sqrt(1 - CPW.k(self) ** 2)
            F = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
        else:
            F = np.log(2 * (1 + np.sqrt(CPW.k(self))) / (1 - np.sqrt(CPW.k(self)))) / np.pi

        return 30 * np.pi / np.sqrt((self.epsilon_r + 1) / 2) / F

    """inductance per unit length"""

    def calc_inductance_geo(self):

        if CPW.k(self) <= 1 / np.sqrt(2):
            kprime = np.sqrt(1 - CPW.k(self) ** 2)
            F = np.pi / np.log(2 * (1 + np.sqrt(kprime)) / (1 - np.sqrt(kprime)))
        else:
            F = np.log(2 * (1 + np.sqrt(CPW.k(self))) / (1 - np.sqrt(CPW.k(self)))) / np.pi

        return pyc.mu_0 / 4 / F

    def calc_inductance_kin(self):  ## zero order approximation
        return self.L_kin_sq / self.width

    def inductance_total(self):
        return CPW.calc_inductance_geo(self) + CPW.calc_inductance_kin(self)

    def calc_inductance_fraction(self):
        return CPW.calc_inductance_geo(self) / (CPW.calc_inductance_geo(self) + CPW.calc_inductance_kin(self))


def _example():

    cpwline = CPW(width=5e-6, spacing=20e-6, epsilon_r=7.5, L_kin_sq=5e-12)
    print('k', cpwline.k())
    print('total capacitance', cpwline.calc_capacitance_total())
    print('characteristic impedance', cpwline.calc_impedance_characteristic())

    print('k', cpwline.k())
    print('total capacitance', cpwline.calc_capacitance_total())
    print('characteristic impedance', cpwline.calc_impedance_characteristic())

    print('geometric inductance', cpwline.calc_inductance_geo())
    print('kinetic inductance', cpwline.calc_inductance_kin())
    print('inductance fraction', cpwline.calc_inductance_fraction())


if __name__ == '__main__':
    _example()
