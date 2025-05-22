# temp for tequila kupccgsd

from tequila.quantumchemistry import Molecule

geomstring="Li 0.0 0.0 0.0\nH 0.0 0.0 1.6"
mol = Molecule(geometry=geomstring)

H = mol.make_hamiltonian()
U = mol.make_ansatz(name="SPA") # or e.g. UpCCGSD

print("Hamiltonian:")
print(H)
print("Ansatz:")
print(U)