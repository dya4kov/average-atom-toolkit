import numpy as np

from pyaatk.thomas_fermi.eos import FreeEnergy as Ftf

from mendeleev import element

kB       = 1.38064853e-23  # J/K
Avogadro = 6.022140857e+23 # N/mol
hbar     = 1.054571817e-34
me       = 9.10938356e-31
echarge  = 1.60217662e-19
ab       = 5.29177210904e-11
ry       = 13.6056930098
aVol     = ab**3           # m^3
eV       = echarge/kB      # K

elem 	 = element('Al')
rho0     = 1000*elem.density   #2700. kg/m^3
rhomin   = 0.01*rho0        # kg/m^3
rhomax   = 1.0*rho0        # kg/m^3
mass     = 1e-3*elem.atomic_weight # 27.e-3 kg/mol
hartree  = 2*ry            # eV
Z        = elem.atomic_number #13.0
V 		 = mass/(Avogadro*rho0 *aVol)
E0       = 0.76874512422 * Z ** (7.0/3.0)

F = Ftf()
F.setZ(Z)

hartree = 13.605693009*2 # eV

# V = np.array([0.1, 1.0, 10.0])
# print("V  = ", V , E0)
Tmin     	= 10.0 # eV
Tmax     	= 1000.0 # eV
NpointsT 	= 100
T 			= np.linspace(Tmin ,Tmax,NpointsT) / hartree# np.array([1.0, 10.0])

#T  = np.array([10])/ hartree

# V = 10 
# T = np.array([1])
# F.setZ(1)
print("rho = ",rho0)
print("T            E_inner_Tf")
file = open("inner_energy_tf.txt", "w")

for T_cur in T: 
	out = "%12.6e %12.6e" % (round(T_cur * hartree), F(V,T_cur) - T_cur * F.DT(V,T_cur))
	file.write(out + "\n")
	print(out)

# print(F(V, T))
# print(F.DV(V, T))
# print(F.DT(V, T))
# print(F.D2V(V, T))
# print(F.DVT(V, T))
# print(F.D2T(V, T))
