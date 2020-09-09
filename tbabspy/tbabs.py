from . import tbabscore
import numpy as np

param_default = np.array([# 0: N_H in 10^{22} cm^{-2}
                          1.0,
                          # 1-17: Abundances He --> Ni
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          # 18: H2 Abundance
                          0.2,
                          # 19: rho_dust (g cm^{-3}), 20: amin (um),
                          # 21: amax (um), 22: dust size dist. power law
                          1.0, 0.025, 0.25, 3.5,
                          # 23-40: depletion factors H --> Ni
                          1.0, 1.0, 0.5, 1.0, 0.6, 1.0, 0.25, 0.2, 0.02,
                          0.1, 0.6, 0.5, 1.0, 0.003, 0.03, 0.3, 0.05, 0.04,
                          # 41: redshift
                          0.0])

def tbabs(E, NH):

    param = param_default.copy()
    param[0] = NH
    
    E = np.atleast_1d(E)

    phot = tbabscore.tbnew(E, param)

    return phot

def ztbabs(E, NH, z):

    param = param_default.copy()
    param[0] = NH
    param[41] = z
    
    E = np.atleast_1d(E)

    phot = tbabscore.tbnew(E, param)

    return phot

def tbfeo(E, NH, AO, AFe, z):

    param = param_default.copy()
    param[0] = NH
    param[4] = AO
    param[15] = AFe
    param[41] = z
    
    E = np.atleast_1d(E)

    phot = tbabscore.tbnew(E, param)

    return phot

def tbgas(E, NH, z):

    param = param_default.copy()

    param[19] = 0.0  # No dust
    param[23:41] = 1.0  # No depletion

    param[0] = NH
    param[41] = z

    E = np.atleast_1d(E)

    phot = tbabscore.tbnew(E, param)

    return phot


