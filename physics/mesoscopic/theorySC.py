import numpy as np


class Superconductivity(object):
    def __init__(self, delta, phi):
        self.delta = delta
        self.phi = phi

    """superconducting order parameter"""
    def calc_op_sc(self):
        return self.delta*np.exp(1j*self.phi)


def _example():

    sc = Superconductivity(delta=1, phi=np.pi/2)
    print('delta', sc.delta)
    print('phi', sc.phi)
    op = sc.calc_op_sc()
    print('order parameter', op)


if __name__ == '__main__':
    _example()