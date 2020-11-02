from atom import Potential, RotatePoints, EnergyLevel, ElectronStates

class Atom:
	def __init__(self):
		self.xU   = Potential()
		self.e    = EnergyLevel()
		self.N    = ElectronStates()
		self.__RP = RotatePoints()

	def setV(self, V):
		self.xU.setV(V)
		self.e.setV(V)
		self.N.setV(V)
		self.__RP.setV(V)

	def setT(self, T):
		self.xU.setT(T)
		self.e.setT(T)
		self.N.setT(T)
		self.__RP.setT(T)

	def setZ(self, Z):
		self.xU.setZ(Z)
		self.e.setZ(Z)
		self.N.setZ(Z)
		self.__RP.setZ(Z)

	def setVTZ(self, V, T, Z):
		self.xU.setVTZ(V, T, Z)
		self.e.setVTZ(V, T, Z)
		self.N.setVTZ(V, T, Z)
		self.__RP.setVTZ(V, T, Z)

	def setTolerance(self, eps):
		self.xU.setTolerance(eps)
		self.e.setTolerance(eps)
		self.N.setTolerance(eps)
		self.__RP.setTolerance(eps)