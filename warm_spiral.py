import numpy as np
import pickle
from galpy import potential
from galpy.df import dehnendf, evolveddiskdf
from scipy.optimize import minimize
from astropy import units as u
import warnings
warnings.filterwarnings('ignore')

def pickle_data(filename, data):
    output = open(filename, 'wb')
    pickle.dump(data, output)
    output.close()

tform = -10 * u.Gyr
tsteady = 3 * u.Gyr
N = 2 # number of arms
alpha = 10 * u.deg  # pitch angle, p in the paper
H = 0.18 * u.kpc  # arbitrary, not specified in the paper
Rs = 0.38

sp = potential.SpiralArmsPotential(amp=1, N=N, alpha=alpha, H=H, Rs=Rs, phi_ref=0)
phi_ref = -minimize(lambda phi: sp.dens(1, 0, phi[0]), x0=0).x[0]

dfc = dehnendf(beta=0, profileParams=(1/3, 1, 40/220), correct=True, niter=20)
lp = potential.LogarithmicHaloPotential(amp=1, normalize=1)

# 30% density perturbation
sp = potential.SpiralArmsPotential(amp=1, N=N, alpha=alpha, H=H, Rs=Rs, phi_ref=phi_ref)
amp30 = abs(lp.dens(R=1, z=0, phi=0) / sp.dens(R=1, z=0, phi=0) * 0.3)

n = 101

omegas = np.linspace(0, 2, 201)
oortA_array = np.empty(len(omegas))
oortB_array = np.empty(len(omegas))
oortC_array = np.empty(len(omegas))
oortK_array = np.empty(len(omegas))

for k in range(len(omegas)):
    sp = potential.SpiralArmsPotential(amp=amp30, N=N, alpha=alpha, H=H, Rs=Rs, omega=omegas[k], phi_ref=phi_ref)
    sp = potential.DehnenSmoothWrapperPotential(pot=sp,tform=tform, tsteady=tsteady)
    edf = evolveddiskdf(dfc, [lp, sp], to=sp._tform)

    oortA_array[k], grid, derivRGrid, derivphiGrid = edf.oortA(R=1, phi=0, gridpoints=n, derivGridpoints=n,
                                 grid=True, derivphiGrid=True, derivRGrid=True, nsigma=4., returnGrids=True)
    oortB_array[k] = edf.oortB(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=4.)
    oortC_array[k] = edf.oortC(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=4.)
    oortK_array[k] = edf.oortK(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=4.)

pickle_data('./pickles/warm_spiralA.pkl', oortA_array)
pickle_data('./pickles/warm_spiralB.pkl', oortB_array)
pickle_data('./pickles/warm_spiralC.pkl', oortC_array)
pickle_data('./pickles/warm_spiralK.pkl', oortK_array)