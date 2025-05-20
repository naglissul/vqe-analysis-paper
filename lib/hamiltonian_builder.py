from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit import QuantumCircuit
from datetime import datetime
from qiskit.circuit import Parameter


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
