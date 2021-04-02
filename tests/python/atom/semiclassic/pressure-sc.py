import math
import numpy as np
import sys
from multiprocessing import Pool
from time import sleep

from mendeleev import element

from pyaatk.atom import SemiclassicAtom as Atom

# from pyaatk.semiclassic import Atom
# from pyaatk.transport.potential import Tabular
# from pyaatk.transport import CrossSection, Kfunction

kB       = 1.38064853e-23  # J/K
Avogadro = 6.022140857e+23 # N/mol
hbar     = 1.054571817e-34
me       = 9.10938356e-31
echarge  = 1.60217662e-19
ab       = 5.29177210904e-11
ry       = 13.6056930098
aVol     = ab**3           # m^3
eV       = echarge/kB      # K

elem     = element('Al')
rho0     = 1e-3*1000*elem.density   #2700. kg/m^3
rhomin   = 0.01*rho0        # kg/m^3
rhomax   = 1.0*rho0        # kg/m^3
mass     = 1e-3*elem.atomic_weight # 27.e-3 kg/mol
hartree  = 2*ry            # eV
Z        = elem.atomic_number #13.0
print('V = ',mass/(Avogadro*rho0*aVol))

Tmin     	 = 1 #90.0 # eV
Tmax     	 = 10**3 # 150.0 # eV
NpointsV 	 = 2
NpointsT 	 = 200
nmax     	 = 15
sigma_energy = 0.0 / hartree


rho = rho0 #*10**np.linspace(-1., 0., NpointsV)
# T = np.linspace(Tmin, Tmax, NpointsT)/hartree # linear scale
T = 10 ** np.linspace(np.log10(Tmin) , np.log10(Tmax), NpointsT)/ hartree # log scale 


[rrho, TT] = np.meshgrid(rho, T)



def Pressure(V, T):
	atom = Atom(V=V, T=T, Z=Z, nmax=nmax,meshSize = 1000, useContinuous = True)#, sigma = sigma_energy
	E0 = 0.76874512422 * Z ** (7.0/3.0)
	r0          = (3.0*V/4.0/math.pi)**(1.0/3.0)
	xmax 		= 1.0
	xmin 		= 1e-3
	x    		= np.linspace(xmin, xmax, 801)**2
	Niterations = 300
	niter 		= 0
	tol 		= 2e-5
	Enew 		= 1
	Eold 		= 0.
	check 		= abs(Eold - Enew)/abs(Eold + Enew)
	P_tf 		= 0.0

	while check > tol and niter < Niterations:
		atom.update(mixing=0.75)
		Eold 	= Enew
		Enew   	= atom.energyFull()
		check 	= abs(Eold - Enew)/abs(Eold + Enew)
		niter 	+= 1
		
	out = "%12.6e %12.6e %12.6e" % (round(T * hartree, 2), niter, atom.discreteLevelsNumber)
	print(out)#Z_calculated
	P_sc = atom.pressure()#

	return [P_tf, P_sc]


## Single process##
# print("rho = ",rho0)
# print("T            E_inner_Tf            E_inner_sc          ")
# file = open("inner_energy_sc.txt", "w")

# for T_cur in T: #np.array([1.0,5.0 ])
# 	V = mass/(Avogadro*rho*aVol)
# 	[E_inner_Tf,E_inner_sc] = Inner_energy(V,T_cur)
# 	out = "%12.6e %12.6e %12.6e" % (round(T_cur * hartree), E_inner_Tf, E_inner_sc)
# 	file.write(out + "\n")
# 	print(out)

## Multi process##
pool = Pool(4)

TT = TT.ravel()
rrho = rrho.ravel()

# print("rho          T            P_Tf         P_sc          ")
print("T            Niter_max    Discrete_levels_number")
sys.stdout.flush()

results = []
for i in range(TT.size):
	rho = rrho[i]
	T = TT[i]
	V = mass/(Avogadro*rho*aVol)
	results.append([i, rho, hartree*T, pool.apply_async(Pressure, args=(V, T))])

ncompleted = 0
ntotal = TT.size
data_items = []
while ncompleted < ntotal:
	sleep(2.0)
	completed = []
	for itask in range(len(results)):
		if results[itask][3].ready():
			number = results[itask][0] 
			rho = results[itask][1]
			T   = results[itask][2]
			[P_Tf, P_sc] = results[itask][3].get()
			out = "%12.6e %12.6e %12.6e %12.6e" % (rho, T, P_Tf, P_sc)
			# print(out)
			data_items.append([number, out])
			sys.stdout.flush()
			completed.append(itask)
			ncompleted += 1
	completed.sort()
	completed.reverse()
	for itask in completed:
		results.pop(itask)

data_items.sort()
data_items = [item[1] for item in data_items]

with open('pressure_atom_sc.txt', 'w') as file:
	for item in data_items:
		file.write("%s\n" % item)
