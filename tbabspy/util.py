from . import tbabscore
import numpy as np

def dgami(s, x):

    return tbabscore.dgami(s, x)


def vernabs(E, Z):

    if np.isscalar(E):
        return tbabscore.vernabs(E, Z)
    else:
        sigma = np.empty(E.shape)
        for i, e in enumerate(E.flat):
            sigma.flat[i] = tbabscore.vernabs(e, Z)
        return sigma


def phfit2(Z, ion, shell, E):
        
    if np.isscalar(E):
        return tbabscore.phfit2(Z, ion, shell, E)
    else:
        sigma = np.empty(E.shape)
        for i, e in enumerate(E.flat):
            sigma.flat[i] = tbabscore.phfit2(Z, ion, shell, e)
        return sigma
