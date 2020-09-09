from . import tbabscore

def tbabs(E, NH):

    params = np.ones(42)
    params[0] = NH
    
    depletWilm = np.array([1.0, 1.0, 0.5, 1.0, 0.6,
                             1.0, 0.25, 0.2, 0.02, 0.1,
                             0.6, 0.5, 1.0, 0.003,
                             0.03, 0.3, 0.05, 0.04])

    param[18] = 0.2    # A_H2
    param[19] = 1.0    # dust density in g/cm^3
    param[20] = 0.025  #min grain size in um
    param[21] = 0.25   #max grain size in um
    param[22] = 3.5    # -PL slope of grain size distribution

    param[23:42] = depletWilm

    params[42] = 0.0

    phot = tbabscore.tbnew(E, params)

    return phot


def dgami(s, x):

    return tbabscore.dgami(s, x)
