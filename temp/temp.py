# temp for tequila kupccgsd

from tequila.quantumchemistry import Molecule

def make_hamiltonian(geomstring:str):
    # geomstring="Li 0.0 0.0 0.0\nH 0.0 0.0 1.6"
    mol = Molecule(geometry=geomstring, verbose=False)

    H = mol.make_hamiltonian()
    return H

def make_upccgsd(geomstring:str, order:int=1):
    # geomstring="Li 0.0 0.0 0.0\nH 0.0 0.0 1.6"
    mol = Molecule(geometry=geomstring, verbose=False)
    U = mol.make_upccgsd_ansatz(order=order)
    return U
