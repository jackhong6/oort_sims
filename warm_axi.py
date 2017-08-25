import pickle
from galpy import potential
from galpy.df import dehnendf, evolveddiskdf
from astropy import units as u
import warnings
warnings.filterwarnings('ignore')

def pickle_data(filename, data):
    output = open(filename, 'wb')
    pickle.dump(data, output)
    output.close()

tform = -10 * u.Gyr
vel_disp = 40 / 220

lp = potential.LogarithmicHaloPotential(amp=1, normalize=1)
dfc = dehnendf(beta=0, profileParams=(1/3, 1, vel_disp), correct=True, niter=20)
edf = evolveddiskdf(dfc, lp, to=tform)

oortA, grid, derivRGrid, derivphiGrid = edf.oortA(R=1, phi=0, gridpoints=101, derivGridpoints=101,
                                 grid=True, derivphiGrid=True, derivRGrid=True, nsigma=6., returnGrids=True)
oortB = edf.oortB(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=6.)
oortC = edf.oortC(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=6.)
oortK = edf.oortK(R=1, phi=0,
                                 grid=grid, derivphiGrid=derivphiGrid, derivRGrid=derivRGrid, nsigma=6.)

print('A={:f} km/s/kpc ; B={:f} km/s/kpc ; C={:f} km/s/kpc ; K={:f} km/s/kpc'.format(
    oortA*220/8, oortB*220/8, oortC*220/8, oortK*220/8))

pickle_data('./pickles/warm_axiA.pkl', oortA)
pickle_data('./pickles/warm_axiB.pkl', oortB)
pickle_data('./pickles/warm_axiC.pkl', oortC)
pickle_data('./pickles/warm_axiK.pkl', oortK)