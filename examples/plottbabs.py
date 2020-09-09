import numpy as np
import matplotlib.pyplot as plt
import tbabspy

Ebins = np.geomspace(0.3, 10, 10000)
E = 0.5*(Ebins[:-1] + Ebins[1:])

abs21 = tbabspy.tbabs(Ebins, 0.1)
abs22 = tbabspy.tbabs(Ebins, 1.0)
abs23 = tbabspy.tbabs(Ebins, 10.0)

fig, ax = plt.subplots(1, 1)
ax.plot(E, abs21, label=r'$N_{\mathrm{H}} = 10^{21}$ cm$^{-2}$')
ax.plot(E, abs22, label=r'$N_{\mathrm{H}} = 10^{22}$ cm$^{-2}$')
ax.plot(E, abs23, label=r'$N_{\mathrm{H}} = 10^{23}$ cm$^{-2}$')

ax.set(xlabel=r'$E$ (keV)', xscale='log', xlim=(Ebins.min(), Ebins.max()),
       ylabel=r'$F_\nu$ obs / $F_\nu$ source', yscale='log',
       ylim=(1.0e-6, 2))
ax.legend()

plt.show()

