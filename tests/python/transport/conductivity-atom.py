import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
from time import sleep

from mendeleev import element

# from pyaatk.semiclassic import Atom
from pyaatk.atom import SemiclassicAtom as Atom

from pyaatk.transport.potential import Tabular
from pyaatk.transport import CrossSection, Kfunction

kB       = 1.38064853e-23  # J/K
Avogadro = 6.022140857e+23 # N/mol
hbar     = 1.054571817e-34
me       = 9.10938356e-31
echarge  = 1.60217662e-19
ab       = 5.29177210904e-11
ry       = 13.6056930098
aVol     = ab**3           # m^3
eV       = echarge/kB      # K


elem = element('Al')
rho0     = 1000*elem.density   #2700. kg/m^3
rhomin   = 0.1*rho0        # kg/m^3
rhomax   = 1.0*rho0        # kg/m^3
mass     = 1e-3*elem.atomic_weight # 27.e-3 kg/mol
hartree  = 2*ry            # eV
Z        = elem.atomic_number #13.0


Tmin     	 = 1.0 # eV
Tmax     	 = 10.0 # eV
NpointsV 	 = 50
NpointsT 	 = 5
nmax     	 = 20
sigma_energy = 1 / hartree
useContinuous= True


# def totalEnergy(atom):
# 	etot = 0.
# 	for n in range(1, nmax):
# 		for l in range(0, n):
# 			enl = atom.energyLevel(n, l)
# 			Nnl = atom.electronStates(n, l)
# 			etot += enl*Nnl
# 	return etot

rho = rho0*10**np.linspace(-2., 0., NpointsV)
T = np.linspace(Tmin, Tmax, NpointsT)/hartree

[rrho, TT] = np.meshgrid(rho, T)

ei   = 1.e-4
ef   = 5.e4
nump = 101
k = 10.0**(0.5*np.linspace(np.log10(ei/ry), np.log10(ef/ry), nump))

def electricConductivity(V, T, ni):
	atom = Atom(V=V, T=T, Z=Z, nmax=nmax, useContinuous = useContinuous)#sigma = sigma_energy

	r0 = (3.0*V/4.0/math.pi)**(1.0/3.0)
	xmax = 1.0
	xmin = 1e-3
	x = np.linspace(xmin, xmax, 601)**2
	Niterations = 50
	niter = 0
	tol = 2e-5
	Enew = 1
	Eold = 0.
	check = abs(Eold - Enew)/abs(Eold + Enew)
	while check > tol and niter < Niterations:
		atom.update(mixing=0.75)
		Eold = Enew
		Enew = atom.energyFull()
		check = abs(Eold - Enew)/abs(Eold + Enew)
		niter += 1

	r = x*r0
	p = -r*atom.U(x)
	scpot = Tabular(r, p, eps=1e-4)
	
	sigma = CrossSection(potential=scpot, rmax=10000., nterms=30)
	csei = sigma(k)*ab*ab # cross-section in SI units
	e  = ry*k*k # energy in eV
	
	ve = np.sqrt(2.*e*echarge/me) # electron velocity m/s
	tau = 1./(ni*ve*csei) # relaxation time in s
	
	taufactor = tau[0] # normalize tau
	taunormalized = tau/taufactor
	
	n32 = 3/2
	n52 = 5/2
	n72 = 7/2

	K32 = Kfunction(n=n32, T=T, M=atom.chemicalPotential)
	K52 = Kfunction(n=n52, T=T, M=atom.chemicalPotential)
	K72 = Kfunction(n=n72, T=T, M=atom.chemicalPotential)
	
	ehartree = e/hartree
	efactor32 = (hartree*echarge)**n32
	efactor52 = (hartree*echarge)**n52
	efactor72 = (hartree*echarge)**n72
	K32 = efactor32*taufactor*K32(ehartree, taunormalized)
	K52 = efactor52*taufactor*K52(ehartree, taunormalized)
	K72 = efactor72*taufactor*K72(ehartree, taunormalized)

	econd = -2.*echarge**2 * math.sqrt(2.*me)/(3.*math.pi**2 * hbar**3)*K32
	tcond = 2. * math.sqrt(2.*me)/(3.*math.pi**2 * hbar**3 * T*hartree*eV)*(-K72 + K52**2/K32)

	return [econd, tcond]

pool = Pool(4)

TT = TT.ravel()
rrho = rrho.ravel()

file = open("results.txt", "w", buffering = 1)

print("rho          T            sigma          kterm          ")
file.write("rho          T            sigma          kterm          \n")
sys.stdout.flush()
# file.flush()
results = []
for i in range(rrho.size):
	rho = rrho[i]
	T = TT[i]
	V = mass/(Avogadro*rho*aVol)
	ni = rho/mass*Avogadro
	results.append([rho, hartree*T, pool.apply_async(electricConductivity, args=(V, T, ni))])


ncompleted = 0
ntotal = rrho.size
while ncompleted < ntotal:
	sleep(1.0)
	completed = []
	for itask in range(len(results)):
		if results[itask][2].ready():
			rho = results[itask][0]
			T   = results[itask][1]
			[econd, tcond] = results[itask][2].get()
			out = "%12.6e %12.6e %12.6e %12.6e" % (rho, T, econd, tcond)
			print(out)
			file.write(out + "\n")
			sys.stdout.flush()
			# file.flush()
			completed.append(itask)
			ncompleted += 1
	completed.sort()
	completed.reverse()
	for itask in completed:
		results.pop(itask)
