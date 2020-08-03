import math
import numpy as np
import matplotlib.pyplot as plt

from pyaatk.TF.atom import Potential as Utf
from pyaatk.TF.eos  import ChemicalPotential as Mtf
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

rho0     = 2700.           # kg/m^3
rhomin   = 0.1*rho0        # kg/m^3
rhomax   = 1.0*rho0        # kg/m^3
mass     = 27.e-3          # kg/mol
hartree  = 2*ry            # eV
Z        = 13.0

Tmin     =  1.0 # eV
Tmax     = 10.0 # eV
NpointsV = 21
NpointsT = 10

rho = rho0*10**np.linspace(-2., 0., NpointsV)
ni = rho/mass*Avogadro # ion concentration N/m^3

V = mass/(Avogadro*rho*aVol)
T = np.linspace(Tmin, Tmax, NpointsT)/hartree

M = Mtf()
M.setZ(Z)
U = Utf()
U.setVTZ(1.0, 1.0, Z)

ei   = 1.e-5
ef   = 1.e3
nump = 101
k = 10.0**(0.5*np.linspace(np.log10(ei/ry), np.log10(ef/ry), nump))

for t in range(NpointsT):
	electricConductivity = []
	U.setT(T[t])
	for v in range(NpointsV):
		U.setV(V[v])
		r0 = (3.0*V[v]/4.0/math.pi)**(1.0/3.0)
		xmax = 1.0
		xmin = 1e-3
		x = np.linspace(xmin, xmax, 400)**2
		r = x*r0
		p = -r*U(x)
		tfpot = Tabular(r, p, eps=1e-4)
	
		sigma = CrossSection(potential=tfpot, rmax=10000., nterms=30)
		csei = sigma(k)*ab*ab # cross-section in SI units
		e  = ry*k*k # energy in eV
		
		ve = np.sqrt(2.*e*echarge/me) # electron velocity m/s
		tau = 1./(ni[v]*ve*csei) # relaxation time in s
	
		taufactor = tau[0] # normalize tau
		taunormalized = tau/taufactor
	
		n = 3/2
		K32 = Kfunction(n=n, T=T[t], M=M(V[v],T[t]))
	
		ehartree = e/hartree
		efactor = (hartree*echarge)**n
		K = efactor*taufactor*K32(ehartree, taunormalized)
		econd = -2.*echarge**2 * math.sqrt(2.*me)/(3.*math.pi**2 * hbar**3)*K
		electricConductivity.append(econd)
		print(v, rho[v], econd)
	
	electricConductivity = np.array(electricConductivity)
	data = np.zeros(electricConductivity.size, dtype=[('rho', 'f8'), ('econd', 'f8')])
	data['rho'] = rho
	data['econd'] = electricConductivity
	np.savetxt('e-conductiviy-TF-Al-' + str(round(hartree*T[t])) + 'eV.txt', data)
	plt.plot(rho, electricConductivity)

plt.xscale('log')
plt.yscale('log')
plt.savefig('e-conductivity-TF-Al.pdf', format='pdf')
