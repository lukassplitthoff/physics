import numpy as np
import scipy.constants as pyc


class SolidState(object):
    @staticmethod
    def energy_zeeman(gl, L, gs, S, B, E_offset):
        return E_offset + pyc.e / 2 / pyc.m_e * (gl * L + gs * S) * B

    @staticmethod
    def energy_zeeman_simple(ms, ge, B, E_offset):
        return E_offset + pyc.e * pyc.hbar / 2 / pyc.m_e * ms * ge * B


def _example():

    e1 = SolidState.energy_zeeman(gl=1, L=0, gs=2, S=1/2, B=np.linspace(0,1, 1001), E_offset=0)
    e2 = SolidState.energy_zeeman_simple(ms=1, ge=2,  B=np.linspace(0, 1, 1001), E_offset=0)
    print('Zeeman energy', e1)
    print('Zeeman energy', e2)


if __name__ == '__main__':
    _example()