import numpy as np


class HybridNW(object):
    def __init__(self, frequency, inductance, resistance, coefficient):
        self.frequency_angular = 2 * np.pi * frequency
        self.inductance = inductance
        self.resistance = resistance
        self.coefficient = coefficient

    def impedance(self):
        impt_inv = 0
        for i in range(len(self.coefficient)):
            impi = self.coefficient[i] * (1j * self.frequency_angular * self.inductance[i] + self.resistance[i])
            impt_inv += 1/impi

        return 1/impt_inv

def _example():

    nw1 = HybridNW(frequency=4e9, inductance=np.array([1e-7]), resistance=np.array([0]), coefficient=np.array([1]))
    print('total impedance NW1', nw1.impedance())
    nw2 = HybridNW(frequency=4e9, inductance=np.array([1e-7, 1.5e-7]), resistance=np.array([0, 1]), coefficient=np.array([1, 0.5]))
    print('total impedance NW2', nw2.impedance())


if __name__ == '__main__':
    _example()