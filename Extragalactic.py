# CRPropa steering file for the propagation from sources to Earth including deflections in random EGMFs
import sys
from crpropa import *
from pylab import *

Run = str(sys.argv[1])
filename_output = 'Extragalactic_'+Run+'_'

# set up random turbulent field
lMax = float(sys.argv[2])/10. * Mpc
Brms = float(sys.argv[3])/10. * nG
alpha = -11./3. #Kolmogorov spectrum

#print(lMax/Mpc)
#print(Brms/nG)
#print(alpha)

mass = int(sys.argv[4])
charge = int(sys.argv[5])

sources = array(['NGC253','NGC4945','Circinus','M83','NGC1808','NGC1068'])
sourceDist = array([3.56,3.72,4.2,4.66,9.08,19.]) * Mpc

sourcePos = Vector3d(0)
gridPoints = 1024
spacing = 60. * kpc
lMin = 2. * spacing
size = (gridPoints - 1.) * spacing
#print('size = '+str(size/Mpc)+' Mpc')
origin = Vector3d(-1., -1., -1.) * size / 2.
minE = 30. * EeV
maxE = 100000. * EeV #1000. * EeV
particles = int(2e3) #2e2 #400000 #4000000
#print(particles)

vgrid = Grid3f(origin, gridPoints, spacing)
initTurbulence(vgrid, Brms, lMin, lMax, alpha)
Bfield = MagneticFieldGrid(vgrid)

#print('Lc = '+str(turbulentCorrelationLength(lMin, lMax, alpha) / Mpc)+' Mpc')  # correlation length

#%%

tolerance = 1.e-4
minStep = 0.1 * kpc
maxStep = 1. * Mpc

# simulation setup
EBL = IRB_Gilmore12
sim = ModuleList()
sim.add(PropagationCK(Bfield, tolerance, minStep, maxStep))
sim.add(PhotoPionProduction(CMB))
sim.add(PhotoPionProduction(EBL))
sim.add(PhotoDisintegration(CMB))
sim.add(PhotoDisintegration(EBL))
sim.add(ElectronPairProduction(CMB))
sim.add(ElectronPairProduction(EBL))
sim.add(NuclearDecay())
sim.add(MinimumEnergy(minE))

# define the observers
for i in range(len(sources)):
    obs = Observer()
    obs.add(ObserverLargeSphere(sourcePos, sourceDist[i]))
    output = HDF5Output(filename_output+sources[i]+'.h5', Output.Everything)
    obs.onDetection(output)
    obs.setDeactivateOnDetection(False)
    #if(i!=len(sources)-1): obs.setDeactivateOnDetection(False)
    #else: obs.setDeactivateOnDetection(True)
    sim.add(obs)

obsBox = Observer()
obsBox.add(ObserverSurface(ParaxialBox(origin,Vector3d(size))))
obsBox.setDeactivateOnDetection(True)
sim.add(obs)

# define the source(s)
source = Source()
source.add(SourcePosition(sourcePos))
source.add(SourceIsotropicEmission())
source.add(SourceParticleType(nucleusId(mass, charge)))
source.add(SourcePowerLawSpectrum(minE, maxE, -1.))

# run simulation
sim.setShowProgress(True)
sim.run(source, particles)
