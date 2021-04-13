import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
from time import sleep

from mendeleev import element

from pyaatk.atom import SemiclassicAtom as Atom


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
rho0     = 0.1*1000*elem.density # 2700. kg/m^3
rhomin   = 0.01*rho0         # kg/m^3
rhomax   = 1.0*rho0          # kg/m^3
mass     = 1e-3*elem.atomic_weight # 27.e-3 kg/mol
hartree  = 2*ry              # eV
Z        = elem.atomic_number#13.0

# print(mass/(Avogadro*rho0*aVol))

Tmin     	 = 10.0 # eV
Tmax     	 = 200.0 # eV
NpointsV 	 = 2
NpointsT 	 = 50#191
nmax     	 = 7#20
sigma_energy = 0.0 / hartree


rho = rho0 #*10**np.linspace(-1., 0., NpointsV)
T = np.linspace(Tmin, Tmax, NpointsT)/hartree

[rrho, TT] = np.meshgrid(rho, T)


def Entropy_Energy(V, T):
	atom = Atom(V=V, T=T, Z=Z, nmax=nmax, meshSize = 1000, useContinuous = True)#, sigma = sigma_energy
	E0 = 0.76874512422 * Z ** (7.0/3.0)

	r0          = (3.0*V/4.0/math.pi)**(1.0/3.0)
	xmax 		= 1.0
	xmin 		= 1e-3
	x    		= np.linspace(xmin, xmax, 801)**2
	Niterations = 50
	niter 		= 0
	tol 		= 2e-5
	Enew 		= 1
	Eold 		= 0.
	check 		= abs(Eold - Enew)/abs(Eold + Enew)

	while check > tol and niter < Niterations:
		atom.update(mixing=0.75)
		Eold 	= Enew
		Enew   	= atom.energyFull()
		# Z_calculated = atom.electronStatesDiscrete() + atom.electronStatesContinuous() 
		check 	= abs(Eold - Enew)/abs(Eold + Enew)
		niter 	+= 1

	# print('T     Niter  Z')
	out = "%12.6e %12.6e %12.6e" % (round(T * hartree, 2), niter, atom.discreteLevelsNumber)
	print(out)#Z_calculated
	Entropy_sc = atom.entropy()
	Inner_energy_sc = atom.internalEnergy()

	return [Entropy_sc, Inner_energy_sc]

pool = Pool(4)

TT = TT.ravel()
rrho = rrho.ravel()

# print("rho          T            S_sc         E_sc          ")
print("T            Niter_max    Discrete_levels_number")
sys.stdout.flush()

results = []
for i in range(TT.size):
	rho = rrho[i]
	T = TT[i]
	V = mass/(Avogadro*rho*aVol)
	results.append([i, rho, hartree*T, pool.apply_async(Entropy_Energy, args=(V, T))])

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
			[S_sc, E_sc] = results[itask][3].get()
			out = "%12.6e %12.6e %12.6e %12.6e" % (rho, T, S_sc, E_sc)
			# print(out)
			data_items.append([number, rho, T, S_sc, E_sc])
			sys.stdout.flush()
			completed.append(itask)
			ncompleted += 1
	completed.sort()
	completed.reverse() # ???
	for itask in completed:
		results.pop(itask)

data_items.sort()

S = [item[3] for item in data_items]
E = [item[4] for item in data_items]
data_items = [ [item[1], item[2] , item[3], item[4] ]  for item in data_items]
# print(data_items)

dT = TT[1] - TT[0]
TdSdT =  TT*np.gradient(S, dT) 
dEdT = np.gradient(E, dT)

plt.plot(TT * hartree, TdSdT, '-x',  label = "T dS / dT")
plt.plot(TT * hartree, dEdT, '-x', label = "dE / dT")
plt.title("Thermodynamic equilibrium")
plt.grid()
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel("T, ev")
# plt.xlim([10, 200])
plt.legend()
plt.show()
print('max(TdSdT - dEdT)= ', max(abs(TdSdT - dEdT)))
print('TdSdT[-1] - dEdT[-1]= ', TdSdT[-1] - dEdT[-1])




# with open('td_consistency.txt', 'w') as file:
# 	for item in data_items:
# 		file.write("%s\n" % item)

