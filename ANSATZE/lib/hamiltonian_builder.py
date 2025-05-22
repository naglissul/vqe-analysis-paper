from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import ParityMapper, JordanWignerMapper, BravyiKitaevMapper
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD



def h_and_uccsd(molecule_string: str, mapper: str) -> tuple[QuantumCircuit, QuantumCircuit]:
    """
    molecule_string: H_2, He, LiH, BeH_2
    mapper: JW, BK, Parity
    """

    driver = PySCFDriver(
        atom=molecule_string,
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()
    if mapper == "JW":
        mapper = JordanWignerMapper()
    elif mapper == "BK":
        mapper = BravyiKitaevMapper()
    elif mapper == "Parity":
        mapper = ParityMapper(num_particles=problem.num_particles)
    else:
        raise ValueError(f"Unknown mapper: {mapper}")
    
    ansatz = UCCSD(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper,
        initial_state=HartreeFock(
            problem.num_spatial_orbitals,
            problem.num_particles,
            mapper,
        ),
    )
    return (mapper.map(problem.hamiltonian.second_q_op()), ansatz)


# ++++++++++++++++++JUNK=----------
# ++++++++++++++++
# ++++++++++++++++

def make_jw_and_uccsd_lih(distance: float = 1.6) -> (QuantumCircuit, QuantumCircuit):
    driver = PySCFDriver(
        atom=f'Li 0 0 0; H 0 0 {distance}',
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()
    mapper = JordanWignerMapper()
    ansatz = UCCSD(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper,
        initial_state=HartreeFock(
            problem.num_spatial_orbitals,
            problem.num_particles,
            mapper,
        ),
    )
    return (mapper.map(problem.hamiltonian.second_q_op()), ansatz)

###=================================

def make_hydrogen_hea() -> QuantumCircuit:
    """    
    Create HEA ansatz circuit for H_2 from this paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9979602/

    Returns:
        QuantumCircuit: The ansatz circuit.
    """
    theta1 = Parameter('θ1')
    theta2 = Parameter('θ2')

    ansatz = QuantumCircuit(2)
    ansatz.ry(theta1, 0)
    ansatz.ry(theta2, 1)
    ansatz.cx(0, 1)

    return ansatz

###============HAMILTONIANS==================

def h2_hamiltonian_parity(distance: float = 0.735) -> QuantumCircuit:
    """
    Create a reduced 2-qubit Hamiltonian for a hydrogen molecule with given H-H distance.

    Args:
        distance (float): The distance between the two hydrogen atoms in angstroms.

    Returns:
        QuantumCircuit: The Hamiltonian as a quantum circuit.
    """
    driver = PySCFDriver(
        atom=f'H .0 .0 .0; H .0 .0 {distance}',
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()

    mapper = ParityMapper(num_particles=problem.num_particles)
    return mapper.map(problem.hamiltonian.second_q_op())

###============HAMILTONIANS==================
def hamiltonian_for_molecule_and_mapper(molecule, mapper) -> QuantumCircuit:
   
    driver = PySCFDriver(
        atom=molecule,
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()

    if mapper == "JW":
        mapper = JordanWignerMapper()
    elif mapper == "BK":
        mapper = BravyiKitaevMapper()
    elif mapper == "Parity":
        mapper = ParityMapper(num_particles=problem.num_particles)
    else:
        raise ValueError(f"Unknown mapper: {mapper}")
    return mapper.map(problem.hamiltonian.second_q_op())

